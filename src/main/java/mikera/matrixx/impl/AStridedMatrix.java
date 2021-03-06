package mikera.matrixx.impl;

import mikera.arrayz.impl.IStridedArray;
import mikera.matrixx.AMatrix;
import mikera.matrixx.Matrixx;
import mikera.vectorz.Vectorz;
import mikera.vectorz.impl.AStridedVector;
import mikera.vectorz.util.ErrorMessages;

public abstract class AStridedMatrix extends AArrayMatrix implements IStridedArray {
	private static final long serialVersionUID = -8908577438753599161L;

	protected AStridedMatrix(double[] data, int rows, int cols) {
		super(data, rows, cols);
	}

	public abstract int getArrayOffset();

	public abstract int rowStride();
	
	public abstract int columnStride();	
	
	@Override
	public abstract AStridedMatrix subMatrix(int rowStart, int rows, int colStart, int cols);
	
	@Override
	public AStridedVector getRowView(int i) {
		return Vectorz.wrapStrided(data, getArrayOffset()+i*rowStride(), cols, columnStride());
	}
	
	@Override
	public AStridedVector getColumnView(int i) {
		return Vectorz.wrapStrided(data, getArrayOffset()+i*columnStride(), rows, rowStride());
	}
	
	@Override
	public AStridedVector getBand(int i) {
		int cs=columnStride();
		int rs=rowStride();
		if ((i>cols)||(i<-rows)) throw new IndexOutOfBoundsException(ErrorMessages.invalidBand(this, i));
		return Vectorz.wrapStrided(data, getArrayOffset()+bandStartColumn(i)*cs+bandStartRow(i)*rs, bandLength(i), rs+cs);
	}
	
	@Override
	public abstract void copyRowTo(int row, double[] dest, int destOffset);
	
	@Override
	public abstract void copyColumnTo(int col, double[] dest, int destOffset);
	
	@Override
	public int[] getStrides() {
		return new int[] {rowStride(), columnStride()};
	}
	
	@Override
	public int getStride(int dimension) {
		switch (dimension) {
		case 0: return rowStride();
		case 1: return columnStride();
		default: throw new IllegalArgumentException(ErrorMessages.invalidDimension(this, dimension));
		}
	}
	
	@Override
	public AMatrix getTranspose() {
		return getTransposeView();
	}
	
	@Override
	public AMatrix getTransposeView() {
		return Matrixx.wrapStrided(getArray(),columnCount(),rowCount(),getArrayOffset(),columnStride(),rowStride());
	}
	
	@Override
	public boolean isPackedArray() {
		return (getArrayOffset()==0)&&(columnStride()==1)&&(rowStride()==columnCount())&&(getArray().length==elementCount());
	}
	
	@Override
	public double[] asDoubleArray() {
		if (isPackedArray()) return getArray();
		return null;
	}
}
