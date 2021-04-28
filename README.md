## climanomaly.m
 <b>climanomaly</b> plots two lines (y vs. x and y vs. ref) and visualises
positive and negative anomalies by shading the area between both lines in
two different colors. This is useful for visualising anomalies of a time
series relative to a climatology. The function can further be used to
plot anomalies relative to a constant baseline or two threshold baselines
(positive anomaly above upper threshold, negative anomaly below lower
threshold).

### Syntax

 <b>climanomaly</b>(x,y,ref)
 <b>climanomaly</b>(...,'top',ColorSpec)
 <b>climanomaly</b>(...,'bottom',ColorSpec)
 <b>climanomaly</b>(...,'mainline','LineSpec')
 <b>climanomaly</b>(...,'refline','LineSpec')
 [hlin,href,htop,hbot] = CLIMANOMALY(...)


### Description

<b>climanomaly</b>(x,y,ref) plots a y vs. x (main line) and y vs. ref (reference
line) and shades areas line values above zero; blue fills the area
between zero and any values below zero.
- To shade anomalies relative to a variable reference (e.g. a
  climatology) specify ref as a vector the length of y.
- To shade anomalies relative to a constant baseline, specify a single
  ref value.
- To shade anomalies relative to an upper and a lower threshold, specify
  two ref values (e.g., let ref be [-0.4 0.5] to shade all values less
  than 0.4 or greater than 0.5).

<b>climanomaly</b>(...,'top',ColorSpec) specifies the top color shading, which
can be described by RGB values or any of the Matlab short-hand color
names (e.g., 'r' or 'red').

<b>climanomaly</b>(...,'bottom',ColorSpec) specifies the bottom shading color.

<b>climanomaly</b>(...,'mainline','LineSpec')
<b>climanomaly</b>(...,'refline','LineSpec')
Specifies line types, plot symbols and colors of the reference line.
LineSpec is a string of characters, e.g. 'b--*'. Refer to the 'plot'
documentation for more options. By default, the main line will be plotted
as a solid black line ('k-') and the reference line as a dotted black
line ('k:').

[hlin,href,htop,hbot] = <b>climanomaly</b>(...) returns the graphics handles of
the main line, top, and bottom plots, respectively.


### Author Info

Jake Weis, University of Tasmania, Institute for Marine and Antarctic
Studies (IMAS), April 2021

This function is based on the <b>anomaly</b> function, written by Chad A.
Greene (<a href="matlab:web('https://github.com/chadagreene/CDT')">Climate Data Toolbox</a>).

Subfunction used: <b>intersections</b> by Douglas M. Schwarz.

See also: plot, boundedline, area, patch, and fill.
