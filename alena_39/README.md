# Клиническая аннотация генетических вариантов

Индивидуальный вариант: `student_39.vcf`.

## 1. Исходные данные и структура VCF

Файл [`student_39.vcf`](student_39.vcf) содержит 50 однонуклеотидных вариантов (SNV) человека в координатах GRCh38. Источник указан в заголовке: тестовый файл Ensembl VEP `Ensembl_VEP_example_homo_sapiens_GRCh38.vcf`, release 116.

| Поле | Значение | В `student_39.vcf` |
|---|---|---|
| CHROM | хромосома | `21` или `22` |
| POS | позиция на хромосоме, 1-based, GRCh38 | например `25587758` |
| ID | идентификатор варианта | rsID у всех 50 записей |
| REF | референсный аллель | один нуклеотид |
| ALT | альтернативный аллель | один нуклеотид, мультиаллельных записей нет |
| QUAL | качество вызова варианта | `.` у всех записей |
| FILTER | результат фильтрации вызова | `.` у всех записей |
| INFO | дополнительные аннотации | `.` у всех записей |

Колонок FORMAT и образцов нет, поэтому генотип и зиготность неизвестны.

Статистика рассчитана скриптом [`vcf_stats.py`](vcf_stats.py), результат сохранён в [`vcf_stats.txt`](vcf_stats.txt).

| Показатель | Значение |
|---|---|
| Всего записей | 50, все SNV |
| Хромосома 21 | 12 |
| Хромосома 22 | 38 |
| Транзиции (A↔G, C↔T) | 40 |
| Трансверсии | 10 |
| Ts/Tv | 4,00 |

## 2. Аннотация

Использованы три источника. Полные результаты и настройки сохранены в репозитории.

### Ensembl VEP

- Release 116, сборка GRCh38.p14, все транскрипты Ensembl/GENCODE без `--pick`.
- Запрошены MANE Select, canonical, HGVS, номера экзонов и интронов, APPRIS, TSL, CCDS и UniProt.
- Сохранены предсказания SIFT, PolyPhen-2, AlphaMissense, REVEL, CADD и SpliceAI, а также известные варианты и клинические статусы из кэша VEP.
- Настройки: [`vep/vep_settings.txt`](vep/vep_settings.txt).
- Полный ответ: [`vep/vep_output.json`](vep/vep_output.json).
- Экспорт по парам «вариант — транскрипт»: [`vep/vep_output.tsv`](vep/vep_output.tsv), 1285 строк.

### ClinVar

- Источник: NCBI E-utilities, обновление ClinVar от 01.11.2026.
- Сначала выполняется поиск по позиции GRCh38, затем запись сопоставляется с вариантом по canonical SPDI: позиция, REF и ALT должны совпасть.
- Записи других аллелей в той же позиции не относятся к варианту и сохраняются отдельно.
- Использована агрегированная germline-классификация уровня VCV и её review status.
- Настройки: [`clinvar/clinvar_settings.txt`](clinvar/clinvar_settings.txt).
- Таблица: [`clinvar/clinvar_results.tsv`](clinvar/clinvar_results.tsv).
- Исходные ответы: [`clinvar/clinvar_output.json`](clinvar/clinvar_output.json).

### gnomAD

- Датасет `gnomad_r4` (gnomAD v4, GRCh38), экзомы, геномы и объединённые данные joint.
- Вариант запрашивается по идентификатору `<chr>-<pos>-<REF>-<ALT>`.
- В сводной таблице используется joint AF = AC/AN; дополнительно указаны максимальная частота среди групп afr, amr, eas, nfe, sas и фильтры gnomAD.
- Настройки: [`gnomad/gnomad_settings.txt`](gnomad/gnomad_settings.txt).
- Таблица: [`gnomad/gnomad_results.tsv`](gnomad/gnomad_results.tsv).
- Полные ответы: [`gnomad/gnomad_output.json`](gnomad/gnomad_output.json).

## 3. Сводная таблица

[`annotated_variants.tsv`](annotated_variants.tsv) содержит 50 строк — по одной на каждый исходный вариант.

Правило выбора основного транскрипта:

1. MANE Select; если их несколько, выбирается транскрипт с наиболее тяжёлым последствием.
2. Если MANE Select нет — белок-кодирующий Ensembl canonical.
3. Иначе — транскрипт с наиболее тяжёлым последствием.

Для всех 50 вариантов найден MANE Select. В колонках `vep_most_severe_any` и `most_severe_transcripts` дополнительно указано самое тяжёлое последствие среди всех транскриптов. Поэтому, например, `stop_gained` во второстепенном транскрипте не подменяет последствие в MANE Select.

В таблице есть координаты и REF>ALT, ген, последствие и impact, выбранный транскрипт и правило выбора, HGVS, предсказатели, rsID, частоты и фильтры gnomAD, наличие записи ClinVar, классификация, review status и VCV. Колонка `vep_clin_sig` сохраняет все аллельные статусы из кэша VEP, включая конфликтующие.

### Основные результаты

| Показатель | Значение |
|---|---|
| Последствия в MANE Select | 25 synonymous, 20 missense, 4 intron, 1 splice_region + synonymous |
| Самое тяжёлое последствие среди всех транскриптов | 28 missense, 21 synonymous, 1 stop_gained |
| Более тяжёлое последствие в другом транскрипте | 9 вариантов |
| Найдены в gnomAD | 48 из 50 |
| AF < 0,1% | 26 |
| 0,1% ≤ AF < 1% | 9 |
| 1% ≤ AF < 5% | 4 |
| AF ≥ 5% | 9 |
| Не найдены в gnomAD | 2 |
| gnomAD `discrepant_frequencies` | 10 |
| Есть запись ClinVar для точного аллеля | 15 из 50 |
| Benign / Benign/Likely benign / Likely benign | 5 / 5 / 2 |
| Uncertain significance | 2 |
| Conflicting classifications of pathogenicity | 1 |

## 4. Отчёт

Краткий отчёт с разбором трёх вариантов и заключением: [`report.md`](report.md).

Разобраны:

- `rs113253213` (`LARGE1`): синонимичная замена с конфликтующими классификациями ClinVar.
- `rs114289284` (`JOSD1`): редкий миссенс-вариант без записей ClinVar и gnomAD, но с неблагоприятными предсказаниями.
- `rs115706909` (`SDF2L1`): `stop_gained` только во второстепенном транскрипте, в MANE Select — синонимичная замена.

## 5. Воспроизведение

Нужен Python 3.8+; сторонние библиотеки не требуются. Для повторного получения аннотаций нужен доступ в интернет. Все уже полученные результаты входят в репозиторий, поэтому локальная пересборка таблицы интернета не требует.

```bash
python vcf_stats.py student_39.vcf > vcf_stats.txt
python vep_annotate.py student_39.vcf
python clinvar_annotate.py student_39.vcf
python gnomad_annotate.py student_39.vcf
python build_table.py vep/vep_output.json clinvar/clinvar_results.tsv gnomad/gnomad_results.tsv annotated_variants.tsv
```
