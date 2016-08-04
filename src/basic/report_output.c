


void writeLaTeXHeader(FILE* fp_latex, char* consensus_contigs_filename) {
    fprintf(fp_latex, "\\documentclass[a4paper]{article}\n\n");
    fprintf(fp_latex, "\\usepackage[english]{babel}\n");
    fprintf(fp_latex, "\\usepackage[utf8x]{inputenc}\n");
    fprintf(fp_latex, "\\usepackage{amsmath}\n");
    fprintf(fp_latex, "\\usepackage{graphicx}\n");
    fprintf(fp_latex, "\\usepackage[colorinlistoftodos]{todonotes}\n");
    fprintf(fp_latex, "\\usepackage{caption}\n");
    fprintf(fp_latex, "\\usepackage{subcaption}\n");
    fprintf(fp_latex, "\n\\title{\\large{MetaCortex %s}}\n\n", consensus_contigs_filename);
    fprintf(fp_latex, "\\begin{document}\n");
    fprintf(fp_latex, "\\maketitle\n\n");
    fprintf(fp_latex, "\\begin{figure}[h]\n");
    fprintf(fp_latex, "\\centering\n");
    fprintf(fp_latex, "\\begin{subfigure}[b]{0.45\\textwidth}\n");
    fprintf(fp_latex, "\\includegraphics[width=\\textwidth]{E_degrees.png}\n");
    fprintf(fp_latex, "\\caption{\\label{fig:frog}E nodes}\n");
    fprintf(fp_latex, "\\end{subfigure}\n");
    fprintf(fp_latex, "\\begin{subfigure}[b]{0.45\\textwidth}\n");
    fprintf(fp_latex, "\\includegraphics[width=\\textwidth]{I_degrees.png}\n");
    fprintf(fp_latex, "\\caption{\\label{fig:frog}I nodes }\n");
    fprintf(fp_latex, "\\end{subfigure}\n");
    fprintf(fp_latex, "\\begin{subfigure}[b]{0.45\\textwidth}\n");
    fprintf(fp_latex, "\\includegraphics[width=\\textwidth]{X_degrees.png}\n");
    fprintf(fp_latex, "\\caption{\\label{fig:frog}X nodes}\n");
    fprintf(fp_latex, "\\end{subfigure}\n");
    fprintf(fp_latex, "\\begin{subfigure}[b]{0.45\\textwidth}\n");
    fprintf(fp_latex, "\\includegraphics[width=\\textwidth]{Y_degrees.png}\n");
    fprintf(fp_latex, "\\caption{\\label{fig:frog}Y nodes }\n");
    fprintf(fp_latex, "\\end{subfigure}\n");
    fprintf(fp_latex, "\\caption{Distribution of nodes in graph}\\label{fig:E_and_I}\n");
    fprintf(fp_latex, "\\end{figure}\n\n\n");
    fprintf(fp_latex, "\\begin{figure}[h]\n");
    fprintf(fp_latex, "\\centering\n");
    fprintf(fp_latex, "\\begin{subfigure}[b]{0.45\\textwidth}\n");
    fprintf(fp_latex, "\\includegraphics[width=\\textwidth]{E_heatmap.png}\n");
    fprintf(fp_latex, "\\caption{\\label{fig:frog}E nodes}\n");
    fprintf(fp_latex, "\\end{subfigure}\n");
    fprintf(fp_latex, "\\begin{subfigure}[b]{0.45\\textwidth}\n");
    fprintf(fp_latex, "\\includegraphics[width=\\textwidth]{E_heatmap_noI.png}\n");
    fprintf(fp_latex, "\\caption{\\label{fig:frog}I nodes }\n");
    fprintf(fp_latex, "\\end{subfigure}\n");
    fprintf(fp_latex, "\\caption{Heatmap of probability of nodes within the graph}\\label{fig:E_and_I}\n");
    fprintf(fp_latex, "\\end{figure}\n\n\n");
}
