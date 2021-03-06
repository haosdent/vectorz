package mikera.matrixx.impl;

import mikera.vectorz.AVector;
import mikera.vectorz.impl.ZeroVector;

/**
 * Abstract base class to represent sparse matrices with a single non-zero band
 * 
 * Unlike ADiagonalMatrix, this matrix need not be square, and may have non-zero values on
 * an arbitrary diagonal.
 * 
 * @author Mike
 *
 */
public abstract class ASingleBandMatrix extends ABandedMatrix {
	private static final long serialVersionUID = -213068993524224396L;

	/**
	 * Override to specify which band of the matrix is nonzero
	 * @return
	 */
	public abstract int nonZeroBand();
	
	/**
	 * Override to specify which band of the matrix is nonzero
	 * @return
	 */
	public abstract AVector getNonZeroBand();
	
	@Override
	public boolean isSymmetric() {
		if (rowCount()!=columnCount()) return false;
		if ((nonZeroBand()==0)||getNonZeroBand().isZero()) return true;
		return false;
	}
	
	@Override
	public boolean isIdentity() {
		return isSquare()&&(nonZeroBand()==0)&&(getNonZeroBand().elementsEqual(1.0));
	}
	
	@Override
	public boolean isZero() {
		return getNonZeroBand().isZero();
	}
	
	@Override
	public boolean isSparse() {
		return true;
	}
	
	@Override 
	public long nonZeroCount() {
		return getNonZeroBand().nonZeroCount();
	}
	
	@Override
	public AVector getBand(int band) {
		if (band==nonZeroBand()) {
			return getNonZeroBand();
		} else {
			return ZeroVector.create(bandLength(band));
		}
	}
	
	// TODO: inner product with single band matrix should be v.fast
}
