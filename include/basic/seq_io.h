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

/************************************************************************
 * seq_io.h
 ************************************************************************/

#ifndef seq_io_h
#define seq_io_h

#define MAX_FIELD_SIZE 20

typedef struct{
    void (* header_parser)(Sequence * seq);
    char * (* get_index)(Sequence * seq);
    char * instrument;
    int run_number;
    char * flowcell_id;
    int lane;
    int tile;
    int x_pos;
    int y_pos;
    int read;
    boolean filtered;
    int control_number;
    char * index;
}casava_sequence_header;

struct{
    char * filename;
    FILE * file;
    SequenceArray * sequences_read_1;
    SequenceArray * sequences_read_2;
    long long written;
    FileFormat format;
    int capacity;
    int size;
    boolean concurrent;
    boolean paird;
}SequenceBufferWriter;

void append_sequence(char * filename, Sequence * seq, FileFormat format);

void append_sequence_fh(FILE * filename, Sequence * seq, FileFormat format);

void * new_sequence_header(sequence_header_type header_type);

#endif
