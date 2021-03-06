package mikera.vectorz.impl;

import java.util.Iterator;

import mikera.arrayz.Arrayz;
import mikera.arrayz.INDArray;
import mikera.arrayz.impl.IStridedArray;
import mikera.matrixx.AMatrix;
import mikera.matrixx.Matrixx;
import mikera.matrixx.impl.StridedMatrix;
import mikera.vectorz.Vector;
import mikera.vectorz.util.ErrorMessages;

/**
 * Abstract base class for vectors backed by a double[] array with a constant stride
 * 
 * The double array can be directly accessed for performance purposes
 * 
 * @author Mike
 */
public abstract class AStridedVector extends ASizedVector implements IStridedArray {
	protected final double[] data;
	
	protected AStridedVector(int length, double[] data) {
		super(length);
		this.data=data;
	}

	private static final long serialVersionUID = -7239429584755803950L;

	public final double[] getArray() {
		return data;
	}
	
	public abstract int getArrayOffset();
	public abstract int getStride();
	
	@Override
	public AStridedVector ensureMutable() {
		return clone();
	}
	
	@Override public double dotProduct(double[] data, int offset) {
		double[] array=getArray();
		int thisOffset=getArrayOffset();
		int stride=getStride();
		int length=length();
		double result=0.0;
		for (int i=0; i<length; i++) {
			result+=array[i*stride+thisOffset]*data[i+offset];
		}
		return result;
	}
	
	@Override
	public double elementSum() {
		int len=length();
		double[] array=getArray();
		int offset=getArrayOffset();
		int stride=getStride();
		double result=0.0;
		for (int i=0; i<len; i++) {
			result+=array[offset+i*stride];
		}		
		return result;
	}
	
	@Override
	public double elementMax(){
		int len=length();
		double[] array=getArray();
		int offset=getArrayOffset();
		int stride=getStride();
		double max = -Double.MAX_VALUE;
		for (int i=0; i<len; i++) {
			double d=array[offset+i*stride];
			if (d>max) max=d;
		}
		return max;
	}
	
	@Override
	public double elementMin(){
		int len=length();
		double[] array=getArray();
		int offset=getArrayOffset();
		int stride=getStride();
		double min = Double.MAX_VALUE;
		for (int i=0; i<len; i++) {
			double d=array[offset+i*stride];
			if (d<min) min=d;
		}
		return min;
	}
	
	@Override
	public INDArray broadcast(int... shape) {
		int dims=shape.length;
		if (dims==0) {
			throw new IllegalArgumentException(ErrorMessages.incompatibleBroadcast(this, shape));
		} else if (dims==1) {
			if (shape[0]!=length()) throw new IllegalArgumentException(ErrorMessages.incompatibleBroadcast(this, shape));
			return this;
		} else if (dims==2) {
			int rc=shape[0];
			int cc=shape[1];
			if (cc!=length()) throw new IllegalArgumentException(ErrorMessages.incompatibleBroadcast(this, shape));
			return Matrixx.wrapStrided(getArray(), rc, cc, getArrayOffset(), 0, getStride());
		}
		if (shape[dims-1]!=length()) throw new IllegalArgumentException(ErrorMessages.incompatibleBroadcast(this, shape));
		int[] newStrides=new int[dims];
		newStrides[dims-1]=getStride();
		return Arrayz.wrapStrided(getArray(),getArrayOffset(),shape,newStrides);
	}
	
	@Override
	public INDArray broadcastLike(INDArray target) {
		if (target instanceof AMatrix) {
			return broadcastLike((AMatrix)target);
		}
		return broadcast(target.getShape());
	}
	
	@Override
	public INDArray broadcastLike(AMatrix target) {
		if (length()==target.columnCount()) {
			return StridedMatrix.wrap(getArray(), target.rowCount(), length(), getArrayOffset(), 0, getStride());
		} else {
			throw new IllegalArgumentException(ErrorMessages.incompatibleBroadcast(this, target));
		}
	}
	
	@Override
	public AStridedVector clone() {
		return Vector.create(this);
	}
	
	@Override
	public void add(Vector v) {
		int length=length();
		if(length!=v.length()) throw new IllegalArgumentException(ErrorMessages.incompatibleShapes(this, v));
		for (int i = 0; i < length; i++) {
			addAt(i,v.data[i]);
		}
	}
	
	@Override
	public void add(double[] data, int offset) {
		int stride=getStride();
		double[] tdata=getArray();
		int toffset=getArrayOffset();
		int length=length();
		for (int i = 0; i < length; i++) {
			tdata[toffset+i*stride]+=data[offset+i];
		}
	}
	
	@Override
	public void addToArray(int offset, double[] destData, int destOffset,int length) {
		int stride=getStride();
		double[] tdata=getArray();
		int toffset=getArrayOffset()+offset*stride;
		for (int i = 0; i < length; i++) {
			destData[destOffset+i]+=tdata[toffset+i*stride];
		}
	}
	
	@Override
	public void addToArray(double[] dest, int destOffset, int destStride) {
		int stride=getStride();
		double[] tdata=getArray();
		int toffset=getArrayOffset();
		for (int i = 0; i < length; i++) {
			dest[destOffset+i*destStride]+=tdata[toffset+i*stride];
		}
	}
	
	@Override
	public double[] asDoubleArray() {
		if (isPackedArray()) return getArray();
		return null;
	}

	@Override
	public boolean isPackedArray() {
		return (getStride()==1)&&(getArrayOffset()==0)&&(getArray().length==length());
	}
	
	@Override
	public int[] getStrides() {
		return new int[] {getStride()};
	}
	
	@Override
	public Iterator<Double> iterator() {
		return new StridedElementIterator(getArray(),getArrayOffset(),length(),getStride());
	}
	
	@Override
	public int getStride(int dimension) {
		switch (dimension) {
		case 0: return getStride();
		default: throw new IllegalArgumentException(ErrorMessages.invalidDimension(this, dimension));
		}
	}
	
	@Override
	public void fill(double value) {
		int stride=getStride();
		double[] array=getArray();
		int di=getArrayOffset();
		for (int i=0; i<length; i++) {
			array[di]=value;
			di+=stride;
		}
	}
	
	@Override
	public boolean equalsArray(double[] data, int offset) {
		int stride=getStride();
		double[] array=getArray();
		int di=getArrayOffset();
		for (int i=0; i<length; i++) {
			if (data[offset+i]!=array[di]) return false;
			di+=stride;
		}
		return true;
	}
}
