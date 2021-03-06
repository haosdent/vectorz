package mikera.matrixx.impl;

import java.util.Arrays;
import java.util.Iterator;

import mikera.arrayz.ISparse;
import mikera.matrixx.AMatrix;
import mikera.matrixx.Matrix;
import mikera.matrixx.Matrixx;
import mikera.randomz.Hash;
import mikera.vectorz.AVector;
import mikera.vectorz.Vector;
import mikera.vectorz.Vectorz;
import mikera.vectorz.impl.RepeatedElementIterator;
import mikera.vectorz.util.DoubleArrays;
import mikera.vectorz.util.ErrorMessages;

/**
 * Lightweight immutable zero matrix class
 */
public final class ZeroMatrix extends ARectangularMatrix implements IFastRows, IFastColumns, ISparse {
	private static final long serialVersionUID = 875833013123277805L;

	@Override public 
	boolean isFullyMutable() {
		return false;
	}
	
	private ZeroMatrix(int rows, int columns) {
		super(rows,columns);
	}
	
	public static ZeroMatrix create(int rows, int columns) {
		return new ZeroMatrix(rows,columns);
	}
	
	public static ZeroMatrix createSameShape(AMatrix a) {
		return new ZeroMatrix(a.rowCount(),a.columnCount());
	}
	
	@Override
	public boolean isSparse() {
		return true;
	}
	
	@Override
	public boolean isSquare() {
		return cols==rows;
	}
	
	@Override
	public boolean isMutable() {
		return false;
	}
	
	@Override
	public boolean isSymmetric() {
		return isSquare();
	}
	
	@Override
	public boolean isDiagonal() {
		return isSquare();
	}
	
	@Override
	public boolean isUpperTriangular() {
		return true;
	}
	
	@Override
	public boolean isLowerTriangular() {
		return true;
	}
	
	@Override
	public boolean isBoolean() {
		return true;
	}
	
	@Override
	public int upperBandwidthLimit() {
		return 0;
	}
	
	@Override
	public int lowerBandwidthLimit() {
		return 0;
	}
	
	@Override
	public void multiply(double factor) {
		// no change - should maybe be an exception because immutable?
	}
	
	@Override
	public AVector getRowView(int row) {
		return Vectorz.createZeroVector(cols);
	}
	
	@Override
	public AVector getColumnView(int col) {
		return Vectorz.createZeroVector(rows);
	}
	
	@Override
	public void copyRowTo(int row, double[] dest, int destOffset) {
		Arrays.fill(dest, destOffset,destOffset+columnCount(),0.0);
	}
	
	@Override
	public void copyColumnTo(int col, double[] dest, int destOffset) {
		Arrays.fill(dest, destOffset,destOffset+rowCount(),0.0);
	}
	
	@Override
	public void addToArray(double[] dest, int offset) {
		// do nothing
	}
	
	@Override
	public void getElements(double[] dest, int destOffset) {
		Arrays.fill(dest, destOffset,destOffset+rowCount()*columnCount(),0.0);
	}

	@Override
	public double determinant() {
		if(isSquare()) throw new UnsupportedOperationException(ErrorMessages.squareMatrixRequired(this));
		return 0.0;
	}
	
	@Override
	public double trace() {
		return 0.0;
	}
	
	@Override
	public double calculateElement(int i, AVector v) {
		assert(i>=0);
		assert(i<rows);
		return 0.0;
	}

	@Override
	public double get(int row, int column) {
		if ((row<0)||(row>=rows)||(column<0)||(column>=cols)) {
			throw new IndexOutOfBoundsException(ErrorMessages.invalidIndex(this, row,column));
		}
		return 0.0;
	}

	@Override
	public void set(int row, int column, double value) {
		throw new UnsupportedOperationException(ErrorMessages.immutable(this));
	}
	
	@Override
	public double unsafeGet(int row, int column) {
		return 0.0;
	}

	@Override
	public void unsafeSet(int row, int column, double value) {
		throw new UnsupportedOperationException(ErrorMessages.immutable(this));
	}
	
	@Override
	public AMatrix clone() {
		return Matrixx.newMatrix(rows, cols);
	}
	
	@Override
	public boolean isZero() {
		return true;
	}
	

	@Override
	public double elementSum() {
		return 0.0;
	}
	
	@Override
	public double elementMax(){
		return 0.0;
	}
	
	@Override
	public double elementMin(){
		return 0.0;
	}
	
	@Override
	public long nonZeroCount() {
		return 0;
	}
	
	@Override 
	public int hashCode() {
		return Hash.zeroVectorHash(cols*rows);
	}
	
	@Override
	public void transform(AVector input, AVector output) {
		assert(output.length()==rows);
		output.fill(0.0);
	}
	
	@Override
	public void transform(Vector input, Vector output) {
		assert(output.length()==rows);
		output.fill(0.0);
	}
	
	@Override
	public boolean isInvertible() {
		return false;
	}
	
	@Override
	public AVector asVector() {
		return Vectorz.createZeroVector(cols*rows);
	}
	
	@Override
	public AMatrix innerProduct(AMatrix m) {
		if (columnCount()!=m.rowCount()) throw new IllegalArgumentException(ErrorMessages.incompatibleShapes(this, m));
		return ZeroMatrix.create(rows, m.columnCount());
	}
	
	@Override
	public ZeroMatrix innerProduct(double a) {
		return this;
	}
	
	@Override 
	public void elementMul(AMatrix m) {
		// do nothing, already zero!
	}

	@Override
	public boolean equals(AMatrix m) {
		if (!isSameShape(m)) return false;
		return m.isZero();
	}
	
	@Override
	public boolean equalsArray(double[] data, int offset) {
		return DoubleArrays.isZero(data, offset, rows*cols);
	}

	@Override
	public ZeroMatrix getTranspose() {
		if (cols==rows) return this;
		return ZeroMatrix.create(cols, rows);
	}
	
	@Override
	public Matrix toMatrix() {
		return Matrix.create(rows, cols);
	}
	
	@Override
	public double[] toDoubleArray() {
		return new double[rows*cols];
	}
	
	@Override
	public AMatrix sparseClone() {
		return Matrixx.createSparse(rows, cols);
	}
	
	@Override
	public Matrix toMatrixTranspose() {
		return Matrix.create(cols, rows);
	}
	
	@Override
	public AVector getLeadingDiagonal() {
		return Vectorz.createZeroVector(Math.min(rows, cols));
	}
	
	@Override
	public AVector getBand(int band) {
		return Vectorz.createZeroVector(bandLength(band));
	}
	
	@Override
	public AMatrix subMatrix(int rowStart, int rows, int colStart, int cols) {
		return ZeroMatrix.create(rows, cols);
	}
	
	@Override
	public Iterator<Double> elementIterator() {
		return new RepeatedElementIterator(cols*rows,0.0);
	}
	
	@Override
	public ZeroMatrix exactClone() {
		return new ZeroMatrix(rows,cols);
	}
}
