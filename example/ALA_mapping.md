
Symbolic Resolution of Labeled Measurements
===========================================

Contents
========

* [SymNmrMs](#symnmrms)
	* [Problem setup](#problem-setup)
	* [Solution](#solution)
	* [System Status](#system-status)
	* [Redundant measurements](#redundant-measurements)
	* [Measurable elementary combinations](#measurable-elementary-combinations)
	* [Least Squares](#least-squares)
  
<a name="symnmrms"></a>
# SymNmrMs
  
<a name="problem-setup"></a>
## Problem setup
  
<a name="used-arguments"></a>
### Used arguments


```
'-t' 'example/ALA_mapping.tsv' '-d' 'example/ALA_measurements.tsv'
```  
<a name="available-measurements"></a>
### Available measurements
  

||HSQC_Ca|HSQC_Cb|HACO|JRES_Ha|JRES_Hb|HACO-DIPSY|GC-MS|LC-MS|
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
|000||||e|t||o|g|
|001||k||e|u||p|h|
|100|||r|e|t||o|h|
|101||k|r|e|u||p|i|
|010|a|||f|t|m|p|h|
|011|b|l||f|u|m|q|i|
|110|c||s|f|t|n|p|i|
|111|d|l|s|f|u|n|q|j|
  
<a name="equation-system"></a>
### Equation system
  
i000 + i001 + i100 + i101 + i010 + i011 + i110 + i111 = 1  
a·(i010 + i011 + i110 + i111) - i010 = 0  
b·(i010 + i011 + i110 + i111) - i011 = 0  
c·(i010 + i011 + i110 + i111) - i110 = 0  
d·(i010 + i011 + i110 + i111) - i111 = 0  
-i001 - i101 + k·(i001 + i011 + i101 + i111) = 0  
-i011 - i111 + l·(i001 + i011 + i101 + i111) = 0  
-i100 - i101 + r·(i100 + i101 + i110 + i111) = 0  
-i110 - i111 + s·(i100 + i101 + i110 + i111) = 0  
e - i000 - i001 - i100 - i101 = 0  
f - i010 - i011 - i110 - i111 = 0  
-i000 - i010 - i100 - i110 + t = 0  
-i001 - i011 - i101 - i111 + u = 0  
-i010 - i011 + m·(i010 + i011 + i110 + i111) = 0  
-i110 - i111 + n·(i010 + i011 + i110 + i111) = 0  
-i000 - i100 + o = 0  
-i001 - i010 - i101 - i110 + p = 0  
-i011 - i111 + q = 0  
g - i000 = 0  
h - i001 - i010 - i100 = 0  
i - i011 - i101 - i110 = 0  
-i111 + j = 0  
<a name="solution"></a>
## Solution
  
<a name="isotopomers"></a>
### Isotopomers

|variable|formula|
| :---: | :---: |
|000|g|
|001|g + h - j·m - n·(-f + j) - t|
|100|-g + o|
|101|i + j·(m + 1) + n·(-f + j) - q|
|010|j·m + n·(-f + j) - o + t|
|011|-j + q|
|110|-j·m - n·(-f + j)|
|111|j|
  
<a name="emu"></a>
### EMU

|variable|Exx|xEx|xxE|EEx|ExE|xEE|EEE|
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
|M+0|2·g + h - j - o + q|e|t|2·g + h - j·m - n·(-f + j) - t|g + j·m + n·(-f + j) - o + t|o|g|
|M+1|-g + i + 2·j + o - q|f|u|-g + i + 2·j·m + 2·n·(-f + j) + t|h - j·(2·m + 1) - 2·n·(-f + j) - p + u|p|h|
|M+2||||f·n|i + j·(m + 2) + n·(-f + j) - q|q|i|
|M+3|||||||j|
  
<a name="cumomers"></a>
### Cumomers

|variable|formula|
| :---: | :---: |
|1xx|-g + i + 2·j + o - q|
|x1x|f|
|xx1|u|
|11x|f·n|
|1x1|i + j·(m + 2) + n·(-f + j) - q|
|x11|q|
|111|j|
  
<a name="system-status"></a>
## System Status


The system is *over-determined*.

The reason is that following measurements were not used: *a*, *b*, *c*, *d*, *e*, *f*, *k*, *l*, *r*, *s*.

Number of redundant measurements: 10.  
<a name="redundant-measurements"></a>
## Redundant measurements

|Measurement|Formula|Methods|
| :---: | :---: | :--- |
|a|(j·m + n·(j + o - q - t) - o + t)/(-o + q + t)|HSQC_Ca: *GC-MS, HACO-DIPSY, JRES_Hb, LC-MS*|
|b|(-j + q)/(-o + q + t)|HSQC_Ca: *GC-MS, JRES_Hb, LC-MS*|
|c|(-j·m - n·(j + o - q - t))/(-o + q + t)|HSQC_Ca: *GC-MS, HACO-DIPSY, JRES_Hb, LC-MS*|
|d|j/(-o + q + t)|HSQC_Ca: *GC-MS, JRES_Hb, LC-MS*|
|e|o - q + u|JRES_Ha: *GC-MS, JRES_Hb*|
|f|-o + q + t|JRES_Ha: *GC-MS, JRES_Hb*|
|k|-(q - u)/u|HSQC_Cb: *GC-MS, JRES_Hb*|
|l|q/u|HSQC_Cb: *GC-MS, JRES_Hb*|
|r|(g - i - j·(m + 1) - n·(j + o - q - t) - o + q)/(g - i - 2·j - o + q)|HACO: *GC-MS, HACO-DIPSY, JRES_Hb, LC-MS*|
|s|-n·(-o + q + t)/(g - i - 2·j - o + q)|HACO: *GC-MS, HACO-DIPSY, JRES_Hb, LC-MS*|
  
<a name="measurable-elementary-combinations"></a>
## Measurable elementary combinations

|N°|Combination|Accumulated|Formula|Methods|
| :---: | :---: | :---: | :---: | :---: |
|1|000|000|g|LC-MS|
|2|001|001|g + h - j·m - n·(-f + j) - t|HACO-DIPSY, JRES_Ha, JRES_Hb, LC-MS|
|3|100|100|-g + o|GC-MS, LC-MS|
|4|101|101|i + j·(m + 1) + n·(-f + j) - q|GC-MS, HACO-DIPSY, JRES_Ha, LC-MS|
|5|010|010|j·m + n·(-f + j) - o + t|GC-MS, HACO-DIPSY, JRES_Ha, JRES_Hb, LC-MS|
|6|011|011|-j + q|GC-MS, LC-MS|
|7|110|110|-j·m - n·(-f + j)|HACO-DIPSY, JRES_Ha, LC-MS|
|8|111|111|j|LC-MS|
  
<a name="least-squares"></a>
## Least Squares
  
<a name="problem-summary"></a>
### Problem summary

|Item|Value|Comment|
| :---: | :---: | :---: |
|number of measurements|21|m|
|system rank (number of statistically defined parameters)|7|p|
|degree of freedom|14|dof=m - p|
|χ²(14)|3.179|Σ((estimated-measured)/sd)²<br/>95% confidence interval is [0; 24]|
|χ²(14)/14|0.227|Reduced χ², χ²(dof)/dof: if ≫ 1 then poor fitting, ~ 1 means 'good' fitting and ≪ 1 means over-fitting or overestimating sd|
|p-value of one-tail χ² test, i.e.<br/>P(χ²(14) > 3.179)|0.999|Value close to 0 (e.g. under 0.05) means poor fitting. Value close to 1 can be an evidence for over-fitting or that sd are overestimated. It can be NaN (not a number) for dof=0|
  
<a name="measured-values"></a>
### Measured values
  

|name|sd|value|estimated|(estimated - value)/sd|
| :---: | :---: | :---: | :---: | :---: |
|a|0.02004|0.254|0.249|-0.25475|
|b|0.02004|0.254|0.244|-0.49063|
|c|0.02004|0.248|0.248|-0.0216|
|d|0.02004|0.244|0.26|0.76697|
|k|0.02|0.497|0.496|-0.04218|
|l|0.02|0.503|0.504|0.04218|
|r|0.03|0.502|0.511|0.29765|
|s|0.03|0.498|0.489|-0.29765|
|e|0.03|0.51|0.497|-0.42222|
|f|0.03|0.49|0.503|0.42222|
|t|0.03|0.496|0.497|0.02311|
|u|0.03|0.504|0.503|-0.02311|
|m|0.03|0.496|0.493|-0.11238|
|n|0.03|0.504|0.507|0.11238|
|o|0.01|0.25|0.247|-0.29363|
|p|0.01|0.5|0.499|-0.05006|
|q|0.01|0.25|0.253|0.34369|
|g|0.01|0.117|0.122|0.533|
|h|0.01|0.355|0.358|0.34705|
|i|0.01|0.387|0.389|0.16109|
|j|0.01|0.141|0.131|-1.04114|
  
<a name="redundant-measurements"></a>
### Redundant measurements

|Measurement|Measured Value|Formula Applied|Δ=FA-MV|
| :---: | :---: | :---: | :---: |
|a|0.254|0.276|0.023|
|b|0.254|0.220|-0.034|
|c|0.248|0.219|-0.028|
|d|0.244|0.285|0.040|
|e|0.51|0.504|-0.006|
|f|0.49|0.496|0.006|
|k|0.497|0.504|0.007|
|l|0.503|0.496|-0.007|
|r|0.502|0.548|0.046|
|s|0.498|0.452|-0.046|
  
<a name="isotopomers"></a>
### Isotopomers

|variable|value|sd|
| :---: | :---: | :---: |
|000|0.122|0.00844|
|001|0.108|0.01277|
|100|0.125|0.01008|
|101|0.142|0.01054|
|010|0.125|0.00867|
|011|0.123|0.00673|
|110|0.125|0.00745|
|111|0.131|0.00584|
  
<a name="cumomers"></a>
### Cumomers

|variable|value|sd|
| :---: | :---: | :---: |
|1xx|0.522|0.01784|
|x1x|0.503|0.01021|
|xx1|0.503|0.0122|
|11x|0.255|0.00822|
|1x1|0.272|0.01368|
|x11|0.253|0.00607|
|111|0.131|0.00584|
  
<a name="emu-values"></a>
### EMU values

|variable|Exx|xEx|xxE|EEx|ExE|xEE|EEE|
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
|M+0|0.478|0.497|0.497|0.23|0.247|0.247|0.122|
|M+1|0.522|0.503|0.503|0.514|0.48|0.499|0.358|
|M+2||||0.255|0.272|0.253|0.389|
|M+3|||||||0.131|
  
<a name="emu-standard-deviations"></a>
### EMU standard deviations

|variable|Exx|xEx|xxE|EEx|ExE|xEE|EEE|
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
|M+0|0.01784|0.01021|0.0122|0.01673|0.01214|0.00748|0.00844|
|M+1|0.01784|0.01021|0.0122|0.0163|0.01818|0.00733|0.00903|
|M+2||||0.00822|0.01368|0.00607|0.00844|
|M+3|||||||0.00584|
  
<a name="measurable-combinations"></a>
### Measurable combinations

|Combination|Accumulated|value|sd|
| :---: | :---: | :---: | :---: |
|000|000|0.122|0.00844|
|001|001|0.108|0.01277|
|100|100|0.125|0.01008|
|101|101|0.142|0.01054|
|010|010|0.125|0.00867|
|011|011|0.123|0.00673|
|110|110|0.125|0.00745|
|111|111|0.131|0.00584|
