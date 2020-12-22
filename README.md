# AlleleSorting

zarovnane sakvencie
ocakavany format, alely su oznacene ako <sample_name>-h1, <sample_name>-h2, <sample_name>-h3, <sample_name>-h4




pouzijeme **prvy skript**, dostane sekvencie roztriedene na A1 A2 B1 B2 a LOG file, ktory popisuje co sa s kazdou vzorkou stalo.
ak sme zvolili REMOVE,  vysledne sekvencie obsahuju len tie, ktore presli cez TRESHOLD  --> exon based analysis
ak sme zvolili MASK, sekvencie su makovane N. vysledne fasta subory obsahuju vsetky sekvencie, tie ktore nepresli treshold su maskovane - identicke. oznacenie je vsak zmenene ako A1-B2, co umoznije konkatenaciu napr cez AMAS


konkatenovane exony do genov treba odtestovat, ci splnaju TRESHOLD   **SKRIPT_2**   - daj 2 moznosti, no defaultne vyhod tie co nesplnaju.
vypise LOG

ak vyhadzujeme, dostavame parcialne datasety. Vsetko popisane v LOG ---> gene based alignments



geny
