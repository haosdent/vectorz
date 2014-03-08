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

package mikera.matrixx.algo.linsol;

import mikera.matrixx.Matrix;
import mikera.matrixx.ops.CommonOps;

/**
 * A matrix can be easily inverted by solving a system with an identify matrix.
 * The only disadvantage of this approach is that additional computations are
 * required compared to a specialized solution.
 * 
 * @author Peter Abeles
 */
public class InvertUsingSolve {

  public static void invert(ILinearSolver solver, Matrix A, Matrix A_inv,
      Matrix storage) {

    if (A.rowCount() != A_inv.rowCount()
        || A.columnCount() != A_inv.columnCount()) {
      throw new IllegalArgumentException(
          "A and A_inv must have the same dimensions");
    }

    CommonOps.setIdentity(storage);

    solver.solve(storage, A_inv);
  }

  public static void invert(ILinearSolver solver, Matrix A, Matrix A_inv) {

    if (A.rowCount() != A_inv.rowCount()
        || A.columnCount() != A_inv.columnCount()) {
      throw new IllegalArgumentException(
          "A and A_inv must have the same dimensions");
    }

    CommonOps.setIdentity(A_inv);

    solver.solve(A_inv, A_inv);
  }
}
