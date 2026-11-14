# L04限定診断の文書差分と根拠所在（2026-11-14）

基準: recovery/r14-20261109 @ 42d467f8355b4c1fa40cdf753c63a302ae84628b。GitHubとDriveのREADMEを直接読取り比較した。GitHub版は通常のプレーン表記、Drive版はコード書式・太字・説明の補足があるが、対象 L04_T01 segment 283 angle 6（0始まり）、径方向弦1.996849 mm、面法線レイ1.983288 mm、L04通常部最小値0.587704685/0.586435343 mm、暫定閾値0.5 mm、L04曖昧1区間・全体36区間、正式未合格は一致。Drive版は「分類確定とせず」と明示し、GitHub版は「曖昧のまま保持」と記すが意味は同じ。

JSONは両版でスキーマ・項目構造・記述量が異なる。GitHub版の decision は短文、Drive版は「raw ray provenance と隣接区間の確認待ち」と明記し、sampling_only や source archive の階層情報も追加している。共通の数値・SHA256・対象・未証明には矛盾なし。ただし同一バイトの複製ではなく、正本の統一は今回行わない。

既知の所在のみ確認: GitHub artifacts/penholder_R14/reports/individual_wall_samples_R14.json は路線別集計・サンプリング設定・例示を含むが、対象25本の各起点、方向、交点距離、交差STL三角形IDを含む独立の生レイ台帳ではない。GitHub同reports/L04_wall_diagnostic/ は recovery.json、README、limited JSONのみ。Drive既知R14 reportsフォルダには回収断面PNG、README、JSON、既存検証ZIP等があるが、25本の生レイ台帳や隣接区間の専用診断ファイルは見つからない。前回recovery.jsonも診断ZIPと中間ファイル未回収と記録。ZIP内部の再展開・広域探索は行っていないため、存在しないことの全域証明ではない。

再開条件: 元の25本について起点座標、方向ベクトル、交点距離、命中面・三角形ID、レイの除外/分類理由を持つ原本ログ（例: L04_T01_segment283_angle6_raw_ray_hits.json）、及び segment 282〜284 の同じ基準の隣接断面・分類根拠（例: L04_T01_segments282_284_adjacent_evidence.json）を、保存元・STL SHA256・生成方法・座標系・判定閾値付きで回収すること。これらは必要データの説明用名称であり、実在ファイル名を断定しない。取得後は既存 STL とSHA照合して1箇所のみ再分類を検討する。

判定: L04曖昧1区間、全体曖昧36区間を維持。局所回収値は全表面最小壁厚の保証ではない。R13/STL/3MF・形状は不変。レイ再計算、造形、通液、実機、発注、課金なし。GitHubのPRとActions実行0件を確認。旧チャットの生存状態は未確認。
