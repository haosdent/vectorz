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

package mikera.matrixx.algo.decompose.eig.impl;

import mikera.matrixx.Matrix;
import mikera.matrixx.algo.decompose.eig.IEigen;
import mikera.matrixx.data.Complex64F;
import mikera.matrixx.ops.MatrixFeatures;

/**
 * Checks to see what type of matrix is being decomposed and calls different
 * eigenvalue decomposition algorithms depending on the results. This primarily
 * checks to see if the matrix is symmetric or not.
 * 
 * 
 * @author Peter Abeles
 */
public class SwitchingEigen implements IEigen {
  IEigen symmetricAlg;
  IEigen generalAlg;
  boolean symmetric;
  // should it compute eigenvectors or just eigenvalues?
  boolean computeVectors;
  Matrix A = Matrix.create(1, 1);
  // tolerance used in deciding if a matrix is symmetric or not
  private double tol;

  /**
   * 
   * @param computeVectors
   * @param tol Tolerance for a matrix being symmetric
   */
  public SwitchingEigen(int matrixSize, boolean computeVectors, double tol) {
    //FIXME
    //symmetricAlg = DecompositionFactory.eig(matrixSize, computeVectors, true);
    //generalAlg = DecompositionFactory.eig(matrixSize, computeVectors, false);
    this.computeVectors = computeVectors;
    this.tol = tol;
  }

  public SwitchingEigen(int matrixSize) {
    this(matrixSize, true, 1e-8);
  }

  @Override
  public int getNumberOfEigenvalues() {
    return symmetric ? symmetricAlg.getNumberOfEigenvalues() : generalAlg
        .getNumberOfEigenvalues();
  }

  @Override
  public Complex64F getEigenvalue(int index) {
    return symmetric ? symmetricAlg.getEigenvalue(index) : generalAlg
        .getEigenvalue(index);
  }

  @Override
  public Matrix getEigenVector(int index) {
    if (!computeVectors)
      throw new IllegalArgumentException(
          "Configured to not compute eignevectors");

    return symmetric ? symmetricAlg.getEigenVector(index) : generalAlg
        .getEigenVector(index);
  }

  @Override
  public boolean decompose(Matrix orig) {
    A.reshape(orig.rowCount(), orig.columnCount());

    symmetric = MatrixFeatures.isSymmetric(A, tol);

    return symmetric ? symmetricAlg.decompose(A) : generalAlg.decompose(A);

  }

  @Override
  public boolean inputModified() {
    // since it doesn't know which algorithm will be used until a matrix is
    // provided make a copy
    // of all inputs
    return false;
  }
}
