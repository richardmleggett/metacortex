/************************************************************************
 *
 * This file is part of MetaCortex
 *
 * Authors:
 *     Richard M. Leggett (richard.leggett@earlham.ac.uk) and
 *     Martin Ayling (martin.ayling@earlham.ac.uk)
 *
 * MetaCortex is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MetaCortex is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MetaCortex.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************
 *
 * This file is modified from source that was part of CORTEX. The
 * original license notice for that is given below.
 *
 ************************************************************************
 *
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 *
 * CORTEX project contacts:
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * Development team:
 *       R. Ramirez-Gonzalez (Ricardo.Ramirez-Gonzalez@bbsrc.ac.uk)
 *       R. Leggett (richard@leggettnet.org.uk)
 *
 ************************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************/

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <file_format.h>


    
static char fileFormatStrings[][12] = {"UNSPECIFIED", "FASTA", "FASTQ", "CTX", "ROCHE", "HASH", "CSFASTA",  "KMERS",};

static char sequenceHeaderStrings[][15] = {"UNKNOWN", "CASAVA_1.8"};

static void stringToUpperCase(char *string){
    int i;
    unsigned long  len =  strlen(string);
    for (i = 0; i < len ; i++){
        if (isalpha(string[i])){
            string[i] = toupper(string[i]);
        }

    }
}

char * file_format_to_string(FileFormat ff){
    return fileFormatStrings[ff];
}

//BEWARE! This transforms the format to upper case!
FileFormat string_to_file_format(char * format){
    stringToUpperCase(format);
    int i;
    for(i = 0; i < FILE_FORMAT_LAST; i++)
    {
        if (strcmp(format, fileFormatStrings[i]) == 0) {
            return i;
        }
    }
    return UNSPECIFIED_FORMAT;
}

char * sequence_header_type_to_string(sequence_header_type sht){
    return sequenceHeaderStrings[sht];

}
sequence_header_type string_to_sequence_header_type(char * format){
    return CASAVA_1_8;
}

