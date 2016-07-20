


void writeLaTeXHeader(FILE* fp_latex, char* consensus_contigs_filename) {
    fprintf(fp_latex, "\\documentclass[a4paper,11pt,oneside]{article}");
    fprintf(fp_latex, "\\usepackage{graphicx}");
    fprintf(fp_latex, "\\usepackage{url}");
    fprintf(fp_latex, "\\usepackage{multirow}");
    fprintf(fp_latex, "\\usepackage{rotating}");
    fprintf(fp_latex, "\\usepackage{color}");
    fprintf(fp_latex, "\\usepackage[compact]{titlesec}");
    fprintf(fp_latex, "\\usepackage[portrait,top=1cm, bottom=2cm, left=1cm, right=1cm]{geometry}");
    fprintf(fp_latex, "\\usepackage{float}");
    fprintf(fp_latex, "\\restylefloat{table}");
    fprintf(fp_latex, "\\begin{document}");
    fprintf(fp_latex, "\\renewcommand*{\\familydefault}{\\sfdefault}");
    fprintf(fp_latex, "\\section*{\\large{MetaCortex %s}}", consensus_contigs_filename);
}


void includeGraphicsIfExists(FILE* fp_latex, char* preTex, char* filename, char* postTex) {
   File f = new File(filename);

   if (f.exists()) {
       fprintf(fp_latex, preTex);
       fprintf(fp_latex, filename);
       fprintf(fp_latex, postTex);
   } else {
       pw.print(" ");
   }
}


MAIN BLOCK
-----------
writeLaTeXHeader(fp_latex, consensus_contigs_filename);

fprintf(fp_latex, ("\\begin{figure}[H]");
fprintf(fp_latex, ("\\centering");
includeGraphicsIfExists(fp_latex, "\\includegraphics[height=3.5cm]{", FIGURE_FILENAME , "}");
fprintf(fp_latex, ("\\end{figure}");
