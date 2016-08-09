

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "report_output.h"


void writeLaTeXHeader(FILE* fp_latex, char* consensus_contigs_filename) {
    fprintf(fp_latex, "\\documentclass[a4paper, 12pt, oneside]{article}\n\n");
    fprintf(fp_latex, "\\usepackage[a4paper,lmargin=4cm,rmargin=2cm,tmargin=2cm,bmargin=2cm]{geometry}\n");
    fprintf(fp_latex, "\\usepackage[english]{babel}\n");
    //fprintf(fp_latex, "\\usepackage[utf8x]{inputenc}\n");
    fprintf(fp_latex, "\\usepackage{amsmath}\n");
    fprintf(fp_latex, "\\usepackage{graphicx}\n");
    //fprintf(fp_latex, "\\usepackage[colorinlistoftodos]{todonotes}\n");
    fprintf(fp_latex, "\\usepackage{caption}\n");
    fprintf(fp_latex, "\\usepackage{subcaption}\n");
    fprintf(fp_latex, "\n\\title{\\large{MetaCortex %s}}\n\n", consensus_contigs_filename);
    fprintf(fp_latex, "\\begin{document}\n");
    fprintf(fp_latex, "\\maketitle\n\n");
    fprintf(fp_latex, "\\begin{figure}[h]\n");
    fprintf(fp_latex, "\\centering\n");
    fprintf(fp_latex, "\\begin{subfigure}[b]{0.45\\textwidth}\n");
    fprintf(fp_latex, "\\includegraphics[width=\\textwidth]{graphs/E_degrees.png}\n");
    fprintf(fp_latex, "\\caption{\\label{fig:E}E nodes}\n");
    fprintf(fp_latex, "\\end{subfigure}\n");
    fprintf(fp_latex, "\\begin{subfigure}[b]{0.45\\textwidth}\n");
    fprintf(fp_latex, "\\includegraphics[width=\\textwidth]{graphs/I_degrees.png}\n");
    fprintf(fp_latex, "\\caption{\\label{fig:I}I nodes }\n");
    fprintf(fp_latex, "\\end{subfigure}\n");
    fprintf(fp_latex, "\\begin{subfigure}[b]{0.45\\textwidth}\n");
    fprintf(fp_latex, "\\includegraphics[width=\\textwidth]{graphs/X_degrees.png}\n");
    fprintf(fp_latex, "\\caption{\\label{fig:X}X nodes}\n");
    fprintf(fp_latex, "\\end{subfigure}\n");
    fprintf(fp_latex, "\\begin{subfigure}[b]{0.45\\textwidth}\n");
    fprintf(fp_latex, "\\includegraphics[width=\\textwidth]{graphs/Y_degrees.png}\n");
    fprintf(fp_latex, "\\caption{\\label{fig:Y}Y nodes }\n");
    fprintf(fp_latex, "\\end{subfigure}\n");
    fprintf(fp_latex, "\\caption{Distribution of nodes in graph}\\label{fig:E_and_I}\n");
    fprintf(fp_latex, "\\end{figure}\n\n\n");
    fprintf(fp_latex, "\\begin{figure}[h]\n");
    fprintf(fp_latex, "\\centering\n");
    fprintf(fp_latex, "\\begin{subfigure}[b]{0.45\\textwidth}\n");
    fprintf(fp_latex, "\\includegraphics[width=\\textwidth]{graphs/E_heatmap.png}\n");
    fprintf(fp_latex, "\\caption{\\label{fig:Heat}Heatmap of probability of node types}\n");
    fprintf(fp_latex, "\\end{subfigure}\n");
    fprintf(fp_latex, "\\begin{subfigure}[b]{0.45\\textwidth}\n");
    fprintf(fp_latex, "\\includegraphics[width=\\textwidth]{graphs/E_heatmap_noI.png}\n");
    fprintf(fp_latex, "\\caption{\\label{fig:Heat without I}Heatmap of probability of node types (with I nodes removed)}\n");
    fprintf(fp_latex, "\\end{subfigure}\n");
    fprintf(fp_latex, "\\caption{Heatmap of probability of nodes within the graph}\\label{fig:E_and_I}\n");
    fprintf(fp_latex, "\\end{figure}\n\\newpage\n");
}
