/**
 * @file    src/psascan_src/psascan.h
 * @author  Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * @section LICENCE
 *
 * This file is part of pSAscan v0.1.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2014-2015
 *   Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#ifndef __PSASCAN_SRC_PSASCAN_H_INCLUDED
#define __PSASCAN_SRC_PSASCAN_H_INCLUDED

#include <string>

namespace psascan_private {

    void pSAscan(std::string input_filename, std::string output_filename,
                 std::string gap_filename, long ram_use,
                 long max_threads,
                 bool verbose, long gap_buf_size = (1L << 21));

}  // namespace psascan_private


// The main function.
static inline void pSAscan(std::string input_filename, std::string output_filename,
        std::string gap_filename, long ram_use, long max_threads,
        bool verbose) {
    psascan_private::pSAscan(input_filename, output_filename,
                             gap_filename, ram_use, max_threads, verbose);
}

#endif  // __PSASCAN_SRC_PSASCAN_H_INCLUDED
