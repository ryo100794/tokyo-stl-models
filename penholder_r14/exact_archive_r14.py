#!/usr/bin/env python3
"""Lossless, checksum-described ZIP transport for R14 binary outputs and inputs.

The archive preserves source bytes. Archive digests and extracted-file digests
have separate identities; neither numerical nor physical qualification is implied.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import stat
import zipfile

MANIFEST = 'ARCHIVE_CONTENTS_R14.json'
SUMS = 'SHA256SUMS.txt'
GITHUB_REGULAR_FILE_LIMIT = 100 * 1024 * 1024


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def safe_name(name):
    p = PurePosixPath(name)
    if not name or p.is_absolute() or '..' in p.parts or '\\' in name or '\n' in name or '\r' in name:
        raise ValueError(f'Unsafe or ambiguous archive member: {name!r}')
    if p.as_posix() != name or name in {'.', MANIFEST, SUMS}:
        raise ValueError(f'Non-canonical or reserved archive member: {name!r}')
    return p


def json_bytes(record):
    return (json.dumps(record, ensure_ascii=False, indent=2) + '\n').encode('utf-8')


def zip_info(name):
    info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
    info.compress_type = zipfile.ZIP_DEFLATED
    info.create_system = 3
    info.external_attr = (stat.S_IFREG | 0o644) << 16
    return info


def verify_archive(archive):
    with zipfile.ZipFile(archive) as z:
        names = z.namelist()
        if len(names) != len(set(names)) or MANIFEST not in names or SUMS not in names:
            raise ValueError('Archive must have unique members and both checksum manifests')
        manifest = json.loads(z.read(MANIFEST))
        rows = manifest['files']
        expected = {r['path'] for r in rows}
        if len(expected) != len(rows) or set(names) != expected | {MANIFEST, SUMS}:
            raise ValueError('Archive manifest does not cover exactly every payload member')
        sums = ''.join(r['sha256'] + '  ' + r['path'] + '\n' for r in rows).encode('utf-8')
        if z.read(SUMS) != sums:
            raise ValueError('The two internal checksum manifests disagree')
        for row in rows:
            safe_name(row['path'])
            info = z.getinfo(row['path'])
            if info.is_dir() or stat.S_ISLNK(info.external_attr >> 16):
                raise ValueError('Only ordinary file members are allowed')
            h = hashlib.sha256(); size = 0
            with z.open(info) as stream:
                for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b''):
                    h.update(chunk); size += len(chunk)
            if size != row['bytes'] or h.hexdigest() != row['sha256']:
                raise ValueError(f'Archive payload differs: {row["path"]}')
        return manifest


def pack(args):
    root, archive, receipt = args.root.resolve(), args.archive.resolve(), args.receipt.resolve()
    if archive.exists() or receipt.exists():
        raise ValueError('Archive and receipt must be new paths')
    selected = {}
    for relative in args.include:
        rel = safe_name(relative)
        path = root.joinpath(*rel.parts)
        if path.is_symlink() or not path.exists() or not path.resolve().is_relative_to(root):
            raise ValueError(f'Invalid payload input: {relative}')
        paths = sorted(path.rglob('*')) if path.is_dir() else [path]
        for candidate in paths:
            if candidate.is_symlink():
                raise ValueError(f'Symlinks are not supported: {candidate}')
            if not candidate.is_file():
                continue
            name = candidate.relative_to(root).as_posix(); safe_name(name)
            if candidate.resolve() in {archive, receipt}:
                raise ValueError('Archive or external receipt cannot be its own payload')
            selected[name] = candidate
    if not selected:
        raise ValueError('No payload files selected')
    rows = [{'path': name, 'bytes': path.stat().st_size, 'sha256': digest(path)}
            for name, path in sorted(selected.items())]
    manifest = {'schema': 'R14-exact-archive-v1', 'role': args.role,
                'payload_bytes_preserved': True, 'files': rows,
                'numerical_geometry_qualification': 'see SHA-bound validation reports',
                'physical_tests_performed': False}
    archive.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as z:
        for row in rows:
            with selected[row['path']].open('rb') as source, z.open(zip_info(row['path']), 'w', force_zip64=True) as target:
                shutil.copyfileobj(source, target, length=4 * 1024 * 1024)
        z.writestr(zip_info(MANIFEST), json_bytes(manifest))
        z.writestr(zip_info(SUMS), ''.join(r['sha256'] + '  ' + r['path'] + '\n' for r in rows).encode('utf-8'))
    checked = verify_archive(archive)
    result = {'schema': 'R14-exact-archive-receipt-v1', 'role': args.role,
              'archive_path': str(archive), 'archive_bytes': archive.stat().st_size,
              'archive_sha256': digest(archive), 'member_count': len(rows),
              'uncompressed_payload_bytes': sum(r['bytes'] for r in rows),
              'all_member_sha256_verified_from_archive': checked == manifest,
              'below_github_regular_file_limit_100_MiB': archive.stat().st_size < GITHUB_REGULAR_FILE_LIMIT,
              'files': rows, 'external_service_upload_performed': False}
    receipt.parent.mkdir(parents=True, exist_ok=True); receipt.write_bytes(json_bytes(result))
    print(json.dumps({k: v for k, v in result.items() if k != 'files'}, ensure_ascii=False))


def extract(args):
    archive, out = args.archive.resolve(), args.out.resolve()
    if args.expected_archive_sha256 and digest(archive) != args.expected_archive_sha256.lower():
        raise ValueError('Archive SHA256 does not match the independently published receipt')
    if out.exists() and any(out.iterdir()):
        raise ValueError('Extraction destination must be new or empty')
    manifest = verify_archive(archive)
    out.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as z:
        for row in manifest['files']:
            relative = safe_name(row['path']); target = out.joinpath(*relative.parts)
            if not target.resolve().is_relative_to(out):
                raise ValueError('Extraction path escapes the destination')
            target.parent.mkdir(parents=True, exist_ok=True)
            with z.open(row['path']) as source, target.open('xb') as stream:
                shutil.copyfileobj(source, stream, length=4 * 1024 * 1024)
            if target.stat().st_size != row['bytes'] or digest(target) != row['sha256']:
                raise ValueError(f'Extracted bytes differ: {row["path"]}')
        (out / MANIFEST).write_bytes(z.read(MANIFEST)); (out / SUMS).write_bytes(z.read(SUMS))
    print(json.dumps({'archive_sha256': digest(archive), 'out': str(out),
                      'payload_files_verified': len(manifest['files']),
                      'all_extracted_bytes_match': True}, ensure_ascii=False))



def split(args):
    archive, out = args.archive.resolve(), args.out.resolve()
    if args.part_bytes <= 0 or args.part_bytes >= GITHUB_REGULAR_FILE_LIMIT:
        raise ValueError('Part size must be positive and smaller than 100 MiB')
    if out.exists() and any(out.iterdir()):
        raise ValueError('Split output must be new or empty')
    verify_archive(archive)
    out.mkdir(parents=True, exist_ok=True); rows = []
    with archive.open('rb') as stream:
        index = 1
        while True:
            payload = stream.read(args.part_bytes)
            if not payload:
                break
            name = archive.name + f'.part{index:03d}'
            target = out / name; target.write_bytes(payload)
            rows.append({'path': name, 'bytes': len(payload), 'sha256': digest(target)})
            index += 1
    manifest = {'schema': 'R14-archive-split-v1', 'archive_name': archive.name,
                'archive_bytes': archive.stat().st_size, 'archive_sha256': digest(archive),
                'assembly': 'binary concatenation in listed order, then exact_archive_r14 extract',
                'parts': rows, 'external_service_upload_performed': False}
    path = out / 'ARCHIVE_PARTS_R14.json'; path.write_bytes(json_bytes(manifest))
    print(json.dumps({'manifest': str(path), 'parts': len(rows), 'archive_sha256': manifest['archive_sha256']}))


def join(args):
    manifest_path, target = args.manifest.resolve(), args.out.resolve()
    manifest = json.loads(manifest_path.read_text())
    if target.exists():
        raise ValueError('Joined archive output must be a new file')
    if len({row['path'] for row in manifest['parts']}) != len(manifest['parts']):
        raise ValueError('Archive part manifest contains duplicate files')
    for row in manifest['parts']:
        name = safe_name(row['path']); path = manifest_path.parent.joinpath(*name.parts)
        if not path.resolve().is_relative_to(manifest_path.parent) or path.is_symlink():
            raise ValueError('Archive part path escapes its directory')
        if path.stat().st_size != row['bytes'] or digest(path) != row['sha256']:
            raise ValueError(f'Archive part checksum mismatch: {row["path"]}')
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open('xb') as stream:
        for row in manifest['parts']:
            with (manifest_path.parent / row['path']).open('rb') as source:
                shutil.copyfileobj(source, stream, length=4 * 1024 * 1024)
    if target.stat().st_size != manifest['archive_bytes'] or digest(target) != manifest['archive_sha256']:
        raise ValueError('Reassembled archive checksum mismatch')
    verify_archive(target)
    print(json.dumps({'archive': str(target), 'archive_sha256': digest(target), 'exact_reassembly_verified': True}))


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='action', required=True)
    q = sub.add_parser('pack'); q.add_argument('--root', type=Path, required=True)
    q.add_argument('--archive', type=Path, required=True); q.add_argument('--receipt', type=Path, required=True)
    q.add_argument('--role', required=True); q.add_argument('--include', nargs='+', required=True)
    q.set_defaults(func=pack)
    q = sub.add_parser('extract'); q.add_argument('--archive', type=Path, required=True)
    q.add_argument('--out', type=Path, required=True); q.add_argument('--expected-archive-sha256')
    q.set_defaults(func=extract)
    q = sub.add_parser('split'); q.add_argument('--archive', type=Path, required=True)
    q.add_argument('--out', type=Path, required=True); q.add_argument('--part-bytes', type=int, default=3 * 1024 * 1024)
    q.set_defaults(func=split)
    q = sub.add_parser('join'); q.add_argument('--manifest', type=Path, required=True)
    q.add_argument('--out', type=Path, required=True); q.set_defaults(func=join)
    a = p.parse_args(); a.func(a)


if __name__ == '__main__':
    main()
