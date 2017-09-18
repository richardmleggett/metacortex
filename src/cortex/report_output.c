

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "report_output.h"
#include "dB_graph.h"
#include "graph_stats.h"
#include "logger.h"


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

// Write a report of various stats out to a latex file
 void writeLaTeXreport(FILE* fp_latex, int MAX_BRANCHES, int COVERAGE_BIN_SIZE, \
  int COVERAGE_BINS, int GRAPH_LOG10_LIMIT, int NUM_BEST_NODES, long int * Contig_Branches, \
  long int * Coverage_Dist, dBGraph * graph, GraphInfo * nodes_in_graph)
{
  int i;
  float percentage;
  char* seq = calloc(256, 1);

  fprintf(fp_latex, "\n\\textbf{Complexity distribution of total graph (X/Y nodes)}\\\\\n");
  for(i=0;i<MAX_BRANCHES;i++){
      fprintf(fp_latex, "%i\\quad %li\\\\\n",i, Contig_Branches[i]);
  }

  fprintf(fp_latex, "\n\\textbf{Coverage dist}\\\\\n");

  // first two lines are for 1, 2-4 cov. after that stick revert to cov bin size
  fprintf(fp_latex, "1\\quad %li\\\\\n", Coverage_Dist[0]);
  fprintf(fp_latex, "2-4\\quad %i\\\\\n", sum_array(Coverage_Dist,1,3));
  if(COVERAGE_BIN_SIZE>4){
      fprintf(fp_latex, ">4\\(\\leq\\)%i\\quad %i\\\\\n", COVERAGE_BIN_SIZE, sum_array(Coverage_Dist,4,COVERAGE_BIN_SIZE-1));
  }
  fprintf(fp_latex, "\\\\\n");

  // now repeat the coverage output, but for every bin
  for(i=0;i<(COVERAGE_BINS*COVERAGE_BIN_SIZE-1);i+=COVERAGE_BIN_SIZE){
      if(COVERAGE_BIN_SIZE>1){
          fprintf(fp_latex, ">%i\\(\\leq\\)%i\\quad %i\\\\\n",(i)*COVERAGE_BIN_SIZE, (i + 1)*COVERAGE_BIN_SIZE, sum_array(Coverage_Dist, i*COVERAGE_BIN_SIZE, (i + 1)*COVERAGE_BIN_SIZE));
      }
      else{
          fprintf(fp_latex, "%i\\quad %li\\\\\n",i+1, Coverage_Dist[i]);
      }
  }
  fprintf(fp_latex, "\\(\\geq\\)%i   \\quad %li\\\\\n",(COVERAGE_BINS)*COVERAGE_BIN_SIZE, Coverage_Dist[(COVERAGE_BINS*COVERAGE_BIN_SIZE)-1]);

  // kmer figures
  fprintf(fp_latex, "\\\\\n\\textbf{kmers}\\\\\nunique\\quad %lli\\quad\\quad total\\quad %lli",graph->unique_kmers,graph->loaded_kmers);
  percentage=(float)(graph->unique_kmers)/(float)(graph->loaded_kmers);
  fprintf(fp_latex, "\\quad\\quad \\%% of total\\quad %.2f\\\\\n", percentage*100);
  percentage=(float)(nodes_in_graph->largest_subgraph)/(float)(graph->unique_kmers);
  fprintf(fp_latex, "\\\\\n\\textbf{subgraphs}\\\\\nlargest subgraph\\quad %i\\quad\\quad  \\%% of total\\quad %.2f\\\\\n",nodes_in_graph->largest_subgraph, percentage*100);
  fprintf(fp_latex, "num subgraphs\\quad %i\\quad num subgraphs\\(>\\)2k\\quad %i\\\\\n", nodes_in_graph->num_subgraphs, nodes_in_graph->num_subgraphs_2k);
  fprintf(fp_latex, "num subgraphs per E\\textsuperscript{6} kmers\\quad %f\\\\\n", (float)(nodes_in_graph->num_subgraphs)/(float)(graph->unique_kmers));
  fprintf(fp_latex, "branches\\quad %i\\quad\\quad per 1000 nodes\\quad %.2f\\\\\n\\\\\n\\textbf{graph size dist:}\\\\\n", nodes_in_graph->branch_nodes_total, (float)(nodes_in_graph->branch_nodes_total)/1000.0);
  for(i=0;i<GRAPH_LOG10_LIMIT;i++){
      fprintf(fp_latex, "\\(\\leq\\)E\\textsuperscript{%i}\\quad %i\\\\\n",i,nodes_in_graph->subgraph_dist[i]);
  }
  fprintf(fp_latex, "\\\\\n\\textbf{highest coverage kmers}\\\\\n");
  for(i=0;i<NUM_BEST_NODES;i++){
      // seq never re-initialised
      binary_kmer_to_seq(&nodes_in_graph->kmer[i], graph->kmer_size, seq);
      fprintf(fp_latex, "%d\\quad %s\\\\\n", nodes_in_graph->best_coverage[i], seq);
  }
  fprintf(fp_latex, "\\\\\n\\textbf{num simple graphs}\\quad %i\\\\\n", (nodes_in_graph->simple_bubbles));
  fprintf(fp_latex, "\n\\end{document}");
}


// Write a report of various stats out to a latex file
 void writeLaTeXreport_to_log_and_screen(int MAX_BRANCHES, int COVERAGE_BIN_SIZE, \
  int COVERAGE_BINS, int GRAPH_LOG10_LIMIT, int NUM_BEST_NODES, long int * Contig_Branches, \
  long int * Coverage_Dist, dBGraph * graph, GraphInfo * nodes_in_graph)
  {
      int i;
      float percentage;
      char* seq = calloc(256, 1);

      log_and_screen_printf("\nComplexity distribution of total graph (X/Y nodes)\n");
      for(i=0;i<MAX_BRANCHES;i++){
          log_and_screen_printf("%i\t%li\n",i, Contig_Branches[i]);
      }

      log_and_screen_printf("\nCoverage dist\n");

      // first two lines are for 1, 2-4 cov. after that stick revert to cov bin size
      log_and_screen_printf("1\t%li\n", Coverage_Dist[0]);
      log_and_screen_printf("2-4\t%i\n", sum_array(Coverage_Dist,1,3));
      if(COVERAGE_BIN_SIZE>4){
          log_and_screen_printf(">4\t%i\t%i\n", COVERAGE_BIN_SIZE, sum_array(Coverage_Dist,4,COVERAGE_BIN_SIZE-1));
      }
      log_and_screen_printf("\n");

      // now repeat the coverage output, but for every bin
      for(i=0;i<(COVERAGE_BINS*COVERAGE_BIN_SIZE-1);i+=COVERAGE_BIN_SIZE){
          if(COVERAGE_BIN_SIZE>1){
              log_and_screen_printf(">%i\t%i\t%i\n",(i)*COVERAGE_BIN_SIZE, (i + 1)*COVERAGE_BIN_SIZE, sum_array(Coverage_Dist, i*COVERAGE_BIN_SIZE, (i + 1)*COVERAGE_BIN_SIZE));
          }
          else{
              log_and_screen_printf("%i\t%li\n",i+1, Coverage_Dist[i]);
          }
      }
      log_and_screen_printf(">=%i\t%li\n",(COVERAGE_BINS)*COVERAGE_BIN_SIZE, Coverage_Dist[(COVERAGE_BINS*COVERAGE_BIN_SIZE)-1]);

      // kmer figures
      log_and_screen_printf("\nkmers\nunique\t%lli\ttotal\t%lli",graph->unique_kmers,graph->loaded_kmers);
      percentage=(float)(graph->unique_kmers)/(float)(graph->loaded_kmers);
      log_and_screen_printf("\t%% of total\t%.2f\n", percentage*100);
      percentage=(float)(nodes_in_graph->largest_subgraph)/(float)(graph->unique_kmers);
      log_and_screen_printf("\nsubgraphs\nlargest subgraph\t%i\t%% of total\t%.2f\n",nodes_in_graph->largest_subgraph, percentage*100);
      log_and_screen_printf("num subgraphs\t%i\tnum subgraphs >2k\t%i\n", nodes_in_graph->num_subgraphs, nodes_in_graph->num_subgraphs_2k);
      log_and_screen_printf("num subgraphs per E^6 kmers\t%f\n", (float)(nodes_in_graph->num_subgraphs)/(float)(graph->unique_kmers));
      log_and_screen_printf("branches\t%i\tper 1000 nodes\t%.2f\n\ngraph size dist:\n", nodes_in_graph->branch_nodes_total, (float)(nodes_in_graph->branch_nodes_total)/1000.0);
      for(i=0;i<GRAPH_LOG10_LIMIT;i++){
          log_and_screen_printf("<E^%i\t%i\n",i,nodes_in_graph->subgraph_dist[i]);
      }
      log_and_screen_printf("\nhighest coverage kmers\n");
      for(i=0;i<NUM_BEST_NODES;i++){
          // seq never re-initialised
          binary_kmer_to_seq(&nodes_in_graph->kmer[i], graph->kmer_size, seq);
          log_and_screen_printf("%d\t%s\n", nodes_in_graph->best_coverage[i], seq);
      }
      log_and_screen_printf("\nnum simple graphs\t%i\n", (nodes_in_graph->simple_bubbles));
  }

int sum_array(long int * array, int first, int last){
    int sum = 0;
    int i;

    for(i=first; i<=last; i++){
        sum += array[i];
    }
    return sum;
}
