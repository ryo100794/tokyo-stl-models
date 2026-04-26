# plateau-citygml-to-stl

Project PLATEAUの東京都23区CityGMLを使い、指定bbox内の建物ポリゴンをSTLへ変換するCLIツールです。Google Mapsの3Dデータや地図タイルは使いません。データソースはG空間情報センターで公開されているProject PLATEAUの東京都23区CityGMLを前提にしています。

初期版は、建物CityGML内の`gml:posList`をストリーミングで読み、リングを三角形ファンでメッシュ化します。LOD1/LOD2の厳密な分類には依存しません。

別のCodexや作業者へ引き継ぐ場合は、`HANDOFF.md`も参照してください。

## セットアップ

Windows PowerShell:

```powershell
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

Linux/WSL:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## PLATEAUデータのダウンロード方法

1. G空間情報センターでProject PLATEAUの東京都23区データを探します。
2. CityGML形式のZIPファイルのダウンロードURLを取得します。
3. 次の`download`コマンドにURLを指定します。

```bash
python plateau_citygml_to_stl.py download --url "<CityGML ZIP URL>" --out data/tokyo23_citygml.zip
```

## 展開

```bash
python plateau_citygml_to_stl.py unzip --zip data/tokyo23_citygml.zip --out data/tokyo23
```

## 小範囲テスト

最初から23区全体や広範囲を処理せず、小さいbboxと`--max-files`で確認してください。

Windows PowerShell:

```powershell
python plateau_citygml_to_stl.py stl --root data/tokyo23 --out output/test_shinjuku.stl --west 139.695 --east 139.710 --south 35.685 --north 35.695 --scale 0.001 --base-thickness 3.0 --max-files 20
```

Linux/WSL:

```bash
python plateau_citygml_to_stl.py stl --root data/tokyo23 --out output/test_shinjuku.stl --west 139.695 --east 139.710 --south 35.685 --north 35.695 --scale 0.001 --base-thickness 3.0 --max-files 20
```

## 本番範囲

Windows PowerShell:

```powershell
python plateau_citygml_to_stl.py stl `
  --root data/tokyo23 `
  --out output/shinjuku_shibuya_shimokita.stl `
  --west 139.655 `
  --east 139.735 `
  --south 35.640 `
  --north 35.715 `
  --scale 0.0005 `
  --base-thickness 3.0
```

Linux/WSL:

```bash
python plateau_citygml_to_stl.py stl \
  --root data/tokyo23 \
  --out output/shinjuku_shibuya_shimokita.stl \
  --west 139.655 \
  --east 139.735 \
  --south 35.640 \
  --north 35.715 \
  --scale 0.0005 \
  --base-thickness 3.0
```

## 23区の本番例

23区全体を1つのSTLにすると非常に重くなります。まず区や駅周辺などに分割して出力することを推奨します。どうしても広範囲を処理する場合は、十分なメモリと時間を確保してください。

```bash
python plateau_citygml_to_stl.py stl --root data/tokyo23 --out output/tokyo23_area.stl --west 139.56 --east 139.93 --south 35.52 --north 35.82 --scale 0.0005 --base-thickness 3.0
```

## scaleの意味

`--scale`は、実寸メートルからSTL単位への倍率です。

- `--scale 0.001`: 1 mを0.001 STL単位にします。STLをmmとして読むスライサーでは、1 mが1 mm相当になります。
- `--scale 0.0005`: 1 mが0.5 mm相当になります。

座標は、採用された点群のバウンディングボックス中心を原点へ移動してからscaleを掛けます。

## 地面プレート

`--base-thickness`を指定すると、生成済みメッシュ全体を覆う長方形の地面プレートを追加します。値はSTL単位です。例えば、スライサー側でSTL単位をmmとして扱う場合、`--base-thickness 3.0`は厚さ3 mm相当です。

## 3Dプリンタ用の注意点

この初期版は、建物ポリゴンを最短経路でSTL化します。出力STLは水密でない可能性があります。3Dプリンタで使う場合は、Blender、MeshLab、PrusaSlicer、Bambu Studioなどで修復確認してください。

`--repair`を付けると、trimeshによる重複面削除、縮退面削除、未参照頂点削除、法線修正を試みます。ただし、建物底面の完全な閉塞までは保証しません。

## よくあるエラーと対処

`empty STL: no polygons matched the bbox`

- bboxがデータ範囲外の可能性があります。
- `--west --east --south --north`を広げてください。
- `--max-files`が小さすぎて対象タイルに届いていない可能性があります。
- 元データの座標系が想定と違う可能性があります。

処理が非常に遅い、またはメモリが足りない:

- まず`--max-files 20`などで確認してください。
- bboxを小さくしてください。
- 複数エリアに分けてSTL化してください。

ZIPのURLで失敗する:

- G空間情報センターのログインや利用条件確認が必要なURLでは、ブラウザで取得した直接ダウンロードURLを使ってください。
- URLが期限付きの場合は、新しいURLを取得してください。

建物がつぶれる、サイズが大きすぎる:

- `--scale`を調整してください。
- スライサー側の単位解釈も確認してください。

## 地下鉄確認用モデル

東京メトロ9路線と都営地下鉄4路線の統合データ生成、地図マッチ用GeoJSON/CSV出力、粗い地下鉄STL出力にも対応しています。

初期テンプレートを作成します。

```bash
python plateau_citygml_to_stl.py subway-init --out-dir subway_data
```

テンプレートには確認用の概略サンプル駅が入っています。実用時は、国土数値情報N02、GTFS/GTFS-JP、確認済み深度CSVで置き換えてください。

正規化データを作成します。

```bash
python plateau_citygml_to_stl.py subway-normalize \
  --railway-root N02-2025 \
  --gtfs-dir gtfs_tokyo_metro \
  --gtfs-dir gtfs_toei \
  --lines-yaml subway_data/subway_lines.yaml \
  --depth-csv subway_data/subway_depth.csv \
  --aliases-csv subway_data/station_aliases.csv \
  --out-dir normalized_subway \
  --west 139.655 \
  --east 139.735 \
  --south 35.640 \
  --north 35.715 \
  --clip-mode bbox
```

地下鉄だけの粗い確認用STLを作成します。

```bash
python plateau_citygml_to_stl.py subway-rough-stl \
  --normalized-dir normalized_subway \
  --out output/subway_rough_with_shafts.stl \
  --tube-out output/subway_rough_all_lines.stl \
  --shaft-out output/subway_station_shafts.stl \
  --scale 0.0005 \
  --tube-radius 1.0 \
  --shaft-radius 1.5 \
  --tube-sides 8 \
  --shaft-sides 8 \
  --vertical-exaggeration 2.0 \
  --z-offset 0.0 \
  --station-shafts \
  --split-by-line \
  --by-line-dir output/subway_rough_by_line \
  --debug-points output/subway_debug_points.ply
```

PLATEAU建物STLと重ねる場合は、建物STL生成時と同じ`origin_json`を指定してください。`origin_json`がない場合、地下鉄データの範囲から自動originを作るため、地下鉄単体の確認には使えますが、建物STLとは原点が一致しません。

```bash
python plateau_citygml_to_stl.py subway-rough-stl \
  --normalized-dir normalized_subway \
  --origin-json origin.json \
  --out output/subway_rough_with_shafts.stl \
  --scale 0.0005 \
  --station-shafts
```

出力される地図マッチ用データは以下です。

- `normalized_subway/normalized_subway_stations.geojson`
- `normalized_subway/normalized_subway_lines.geojson`
- `normalized_subway/normalized_subway_segments.geojson`
- `normalized_subway/normalized_subway_stations_epsg6677.csv`
- `normalized_subway/normalized_subway_track_points_epsg6677.csv`
- `normalized_subway/subway_data_quality_report.md`

`depth_m`は地表から駅・線路相当位置までの深さです。DEMなしでは`surface_z_m=0.0`、`station_z_m=-depth_m`として扱います。立坑は駅位置確認用で、実際の換気塔、避難立坑、駅出入口の正確な位置を表しません。粗いSTLは地図マッチと可視化用であり、正確なトンネル中心線、連続縦断勾配、施工図面レベルの地下構造物は再現しません。

確認はBlender、MeshLab、PrusaSlicer、Bambu Studioなどで行ってください。
