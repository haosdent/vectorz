/*
 * Copyright (c) 2009-2014, Peter Abeles. All Rights Reserved.
 *
 * This file is part of Efficient Java Matrix Library (EJML).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package mikera.matrixx.algo.misc;

import mikera.matrixx.Matrix;

/**
 * Low level transpose algorithms. No sanity checks are performed. Take a look
 * at BenchmarkTranspose to see which one is faster on your computer.
 * 
 * @author Peter Abeles
 */
public class TransposeAlgs {

  /**
   * In-place transpose for a square matrix. On most architectures it is faster
   * than the standard transpose algorithm, but on most modern computers it's
   * slower than block transpose.
   * 
   * @param mat The matrix that is transposed in-place. Modified.
   */
  public static void square(Matrix mat) {
    int index = 1;
    int indexEnd = mat.columnCount();
    for (int i = 0; i < mat.rowCount(); i++, index += i + 1, indexEnd +=
        mat.columnCount()) {
      int indexOther = (i + 1) * mat.columnCount() + i;
      for (; index < indexEnd; index++, indexOther += mat.columnCount()) {
        double val = mat.get(index);
        mat.set(index, mat.get(indexOther));
        mat.set(indexOther, val);
      }
    }
  }

  /**
   * Performs a transpose across block sub-matrices. Reduces the number of cache
   * misses on larger matrices.
   * 
   * *NOTE* If this is beneficial is highly dependent on the computer it is run
   * on. e.g: - Q6600 Almost twice as fast as standard. - Pentium-M Same speed
   * and some times a bit slower than standard.
   * 
   * @param A Original matrix. Not modified.
   * @param A_tran Transposed matrix. Modified.
   * @param blockLength Length of a block.
   */
  public static void block(Matrix A, Matrix A_tran, final int blockLength) {
    for (int i = 0; i < A.rowCount(); i += blockLength) {
      int blockHeight = Math.min(blockLength, A.rowCount() - i);

      int indexSrc = i * A.columnCount();
      int indexDst = i;

      for (int j = 0; j < A.columnCount(); j += blockLength) {
        int blockWidth = Math.min(blockLength, A.columnCount() - j);

        // int indexSrc = i*A.columnCount() + j;
        // int indexDst = j*A_tran.columnCount() + i;

        int indexSrcEnd = indexSrc + blockWidth;
        // for( int l = 0; l < blockWidth; l++ , indexSrc++ ) {
        for (; indexSrc < indexSrcEnd; indexSrc++) {
          int rowSrc = indexSrc;
          int rowDst = indexDst;
          int end = rowDst + blockHeight;
          // for( int k = 0; k < blockHeight; k++ , rowSrc += A.columnCount() )
          // {
          for (; rowDst < end; rowSrc += A.columnCount()) {
            // faster to write in sequence than to read in sequence
            A_tran.set(rowDst++, A.get(rowSrc));
          }
          indexDst += A_tran.columnCount();
        }
      }
    }
  }

  /**
   * A straight forward transpose. Good for small non-square matrices.
   * 
   * @param A Original matrix. Not modified.
   * @param A_tran Transposed matrix. Modified.
   */
  public static void standard(Matrix A, Matrix A_tran) {
    int index = 0;
    for (int i = 0; i < A_tran.rowCount(); i++) {
      int index2 = i;

      int end = index + A_tran.columnCount();
      while (index < end) {
        A_tran.set(index++, A.get(index2));
        index2 += A.columnCount();
      }
    }
  }
}
