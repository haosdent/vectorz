package mikera.matrixx.impl;

import mikera.matrixx.AMatrix;
import mikera.vectorz.AVector;
import mikera.vectorz.impl.AMatrixViewVector;
import mikera.vectorz.impl.MatrixIndexScalar;
import mikera.vectorz.impl.Vector0;
import mikera.vectorz.util.ErrorMessages;

/**
 * Vector class representing a view of a matrix band
 * @author Mike
 *
 */
@SuppressWarnings("serial")
public final class MatrixBandView extends AMatrixViewVector {
	private int band;
	
	private MatrixBandView(AMatrix source, int band) {
		super(source,source.bandLength(band));
		this.band=band;
	}

	public static AVector create(AMatrix source, int band) {
		int rc=source.rowCount();
		int cc=source.columnCount();
		if ((band>cc)||(band<-rc)) throw new IllegalArgumentException(ErrorMessages.invalidBand(source,band));
		if ((band==cc)||(band==-rc)) return Vector0.INSTANCE;
		return new MatrixBandView(source,band);
	}
	
	@Override
	public void addToArray(double[] data, int offset) {
		int r=source.bandStartRow(band);
		int c=source.bandStartColumn(band);
		for (int i=0; i<length; i++) {
			data[offset+i]+=source.unsafeGet(r+i, c+i);
		}
	}
	
	@Override
	public void getElements(double[] data, int offset) {
		int r=source.bandStartRow(band);
		int c=source.bandStartColumn(band);
		for (int i=0; i<length; i++) {
			data[offset+i]=source.unsafeGet(r+i, c+i);
		}
	}
	
	@Override
	protected int calcRow(int i) {
		return (band<0)?i-band:i;
	}
	@Override
	protected int calcCol(int i) {
		return (band>0)?i+band:i;
	}
	
	@Override
	public MatrixIndexScalar slice(int i) {
		if ((i<0)||(i>=length)) throw new IndexOutOfBoundsException(ErrorMessages.invalidIndex(this, i));
		return MatrixIndexScalar.wrap(source, calcRow(i), calcCol(i));
	}

	@Override
	public AVector exactClone() {
		return new MatrixBandView(source.exactClone(),band);
	}

}
