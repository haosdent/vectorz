package mikera.matrixx.impl;

import mikera.vectorz.AVector;
import mikera.vectorz.Vector;
import mikera.vectorz.Vectorz;
import mikera.vectorz.impl.ZeroVector;
import mikera.vectorz.util.ErrorMessages;
import mikera.vectorz.util.VectorzException;

/**
 * Sparse banded matrix implementation.
 * @author Mike
 *
 */
public class BandedMatrix extends ABandedMatrix {
	private static final long serialVersionUID = -4014060138907872914L;

	private final int minBand;
	private final int maxBand;
	
	private final AVector[] bands;
	private final int rowCount;
	private final int columnCount;

	private BandedMatrix(int rc, int cc, int minBand, AVector[] bands) {
		this.rowCount=rc;
		this.columnCount=cc;
		this.bands=bands;
		this.minBand=minBand;
		this.maxBand=minBand+bands.length-1;
	}
	
	public static BandedMatrix create(int rowCount, int columnCount, int minBand, int maxBand) {
		AVector[] bands=new AVector[maxBand-minBand+1];
		for (int i=minBand; i<=maxBand; i++) {
			bands[i-minBand]=Vector.createLength(bandLength(rowCount,columnCount,i));
		}
		return new BandedMatrix(rowCount,columnCount,minBand,bands);
	}
	
	public static BandedMatrix wrap(int rowCount, int columnCount, int minBand, int maxBand, AVector... bands) {
		if (bands.length!=(maxBand-minBand+1)) throw new IllegalArgumentException("Wrong number of bands: "+bands.length);
		for (int i=minBand; i<=maxBand; i++) {
			AVector b=bands[i-minBand];
			if (b.length()!=bandLength(rowCount,columnCount,i)) {
				throw new IllegalArgumentException("Incorrect length of band "+ i +", was given: "+b.length());
			}
		}
		return new BandedMatrix(rowCount,columnCount,minBand,bands);
	}
	
	@Override
	public int upperBandwidthLimit() {
		return maxBand;
	}

	@Override
	public int lowerBandwidthLimit() {
		return -minBand;
	}
	
	@Override
	public BandedMatrix getTranspose() {
		AVector[] nbands=new AVector[bands.length];
		for (int i=0; i<(-minBand+maxBand+1); i++) {
			nbands[i]=bands[bands.length-1-i];
		}
		return BandedMatrix.wrap(columnCount(), rowCount(), -maxBand, -minBand,nbands);
	}

	@Override
	public AVector getBand(int band) {
		if ((band>=minBand)&&(band<=maxBand)) return bands[band-minBand];
		if ((band>=-rowCount)&&(band<=columnCount)) return Vectorz.createZeroVector(bandLength(band));
		throw new IndexOutOfBoundsException(ErrorMessages.invalidBand(this, band));
	}

	@Override
	public int rowCount() {
		return rowCount;
	}

	@Override
	public int columnCount() {
		return columnCount;
	}
	
	@Override
	public boolean isView() {
		return true;
	}
	
	@Override
	public double get(int row, int column) {
		if ((row<0)||(column<0)||(row>=rowCount)||(column>=columnCount)) throw new IndexOutOfBoundsException("["+row+","+column+"]");
		return getBand(bandIndex(row,column)).unsafeGet(bandPosition(row,column));
	}

	@Override
	public void set(int row, int column, double value) {
		getBand(bandIndex(row,column)).set(bandPosition(row,column),value);
	}
	
	public void addAt(int i, int j, double d) {
		int band=j-i;
		AVector b=getBand(band);
		b.addAt(bandIndex(i,j), d);
	}

	@Override
	public BandedMatrix exactClone() {
		BandedMatrix b=new BandedMatrix(rowCount,columnCount,minBand,bands.clone());
		for (int i=minBand; i<=maxBand; i++) {
			b.bands[i-minBand]=b.bands[i-minBand].exactClone();
		}
		return b;
	}
	
	@Override
	public void transform(AVector source, AVector dest) {
		if (!(dest instanceof Vector)) {
			super.transform(source, dest);
		} else if ((source instanceof Vector )) {
			transform ((Vector)source, (Vector)dest);
		} else {
			Vector t=(Vector)dest;
			t.fill(0.0);
			for (int i=minBand; i<=maxBand; i++) {
				AVector b=getBand(i);
				b.addProductToArray(1.0, 0, source, Math.max(i, 0), t.data, Math.max(-i, 0), bandLength(i));
			}
		}
	}
	
	@Override
	public void transform(Vector source, Vector dest) {
		Vector t=dest;
		t.fill(0.0);
		for (int i=minBand; i<=maxBand; i++) {
			AVector b=getBand(i);
			b.addProductToArray(1.0, 0, source, Math.max(i, 0), t.data, Math.max(-i, 0), bandLength(i));
		}		
	}
	
	@Override public void validate() {
		super.validate();
		if (minBand!=-lowerBandwidthLimit()) throw new VectorzException("Mismatched lower limit: "+minBand);
		if (maxBand!=upperBandwidthLimit()) throw new VectorzException("Mismatched upper limit: "+maxBand);
	}
}
