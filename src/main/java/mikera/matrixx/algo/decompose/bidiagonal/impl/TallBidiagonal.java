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

package mikera.matrixx.algo.decompose.bidiagonal.impl;

import mikera.matrixx.Matrix;
import mikera.matrixx.UtilEjml;
import mikera.matrixx.algo.decompose.bidiagonal.IBidiagonal;
import mikera.matrixx.algo.decompose.qr.IQRP;
import mikera.matrixx.algo.decompose.qr.impl.HouseholderColumnQRP;
import mikera.matrixx.ops.CommonOps;

/**
 * <p>
 * {@link IBidiagonal}
 * specifically designed for tall matrices. First step is to perform QR
 * decomposition on the input matrix. Then R is decomposed using a bidiagonal
 * decomposition. By performing the bidiagonal decomposition on the smaller
 * matrix computations can be saved if m/n > 5/3 and if U is NOT needed.
 * </p>
 * 
 * <p>
 * A = [Q<sub>1</sub> Q<sub>2</sub>][U1 0; 0 I] [B1;0] V<sup>T</sup><br>
 * U=[Q<sub>1</sub>*U1 Q<sub>2</sub>]<br>
 * B=[B1;0]<br>
 * A = U*B*V<sup>T</sup>
 * </p>
 * 
 * <p>
 * A QRP decomposition is used internally. That decomposition relies an a fixed
 * threshold for selecting singular values and is known to be less stable than
 * SVD. There is the potential for a degregation of stability by using
 * BidiagonalDecompositionTall instead of BidiagonalDecomposition. A few simple
 * tests have shown that loss in stability to be insignificant.
 * </p>
 * 
 * <p>
 * See page 404 in "Fundamentals of Matrix Computations", 2nd by David S.
 * Watkins.
 * </p>
 * 
 * 
 * @author Peter Abeles
 */
// TODO optimize this code
public class TallBidiagonal implements IBidiagonal {
  IQRP decompQRP = new HouseholderColumnQRP(); // todo this should be passed in
  IBidiagonal decompBi =
      new RowBidiagonal();
  Matrix B = Matrix.create(1, 1);
  // number of rows
  int m;
  // number of column
  int n;
  // min(m,n)
  int min;

  @Override
  public void getDiagonal(double[] diag, double[] off) {
    diag[0] = B.get(0);
    for (int i = 1; i < n; i++) {
      diag[i] = B.get(i, i);
      off[i - 1] = B.get(i - 1, i);
    }
  }

  @Override
  public Matrix getB(Matrix B, boolean compact) {
    B = RowBidiagonal.handleB(B, compact, m, n, min);

    B.set(0, 0, this.B.get(0, 0));
    for (int i = 1; i < min; i++) {
      B.set(i, i, this.B.get(i, i));
      B.set(i - 1, i, this.B.get(i - 1, i));
    }
    if (n > m)
      B.set(min - 1, min, this.B.get(min - 1, min));

    return B;
  }

  @Override
  public Matrix getU(Matrix U, boolean transpose,
      boolean compact) {
    U = RowBidiagonal.handleU(U, false, compact, m, n, min);

    if (compact) {
      // U = Q*U1
      Matrix Q1 = decompQRP.getQ(true);
      Matrix U1 = decompBi.getU(null, false, true);
      CommonOps.mult(Q1, U1, U);
    } else {
      // U = [Q1*U1 Q2]
      //FIXME
      Matrix Q = decompQRP.getQ(false);
      Matrix U1 = decompBi.getU(null, false, true);
      Matrix Q1 = CommonOps.extract(Q, 0, Q.rowCount(), 0, min);
      Matrix tmp = Matrix.create(Q1.rowCount(), U1.columnCount());
      CommonOps.mult(Q1, U1, tmp);
      CommonOps.insert(tmp, Q, 0, 0);
    }

    if (transpose)
      CommonOps.transpose(U);

    return U;
  }

  @Override
  public Matrix getV(Matrix V, boolean transpose,
      boolean compact) {
    return decompBi.getV(V, transpose, compact);
  }

  @Override
  public boolean decompose(Matrix orig) {

    decompQRP
        .setSingularThreshold(CommonOps.elementMaxAbs(orig) * UtilEjml.EPS);
    if (!decompQRP.decompose(orig)) {
      return false;
    }

    m = orig.rowCount();
    n = orig.columnCount();
    min = Math.min(m, n);
    B.reshape(min, n);

    //FIXME
    decompQRP.getR(true);

    // apply the column pivots.
    // TODO this is horribly inefficient
    Matrix result = Matrix.create(min, n);
    Matrix P = decompQRP.getPivotMatrix(null);
    CommonOps.multTransB(B, P, result);
    B.set(result);

    return decompBi.decompose(B);
  }

  @Override
  public boolean inputModified() {
    return decompQRP.inputModified();
  }
}
