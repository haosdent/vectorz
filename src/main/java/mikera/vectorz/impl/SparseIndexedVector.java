package mikera.vectorz.impl;

import java.util.Arrays;

import mikera.indexz.Index;
import mikera.vectorz.AVector;
import mikera.vectorz.ArrayVector;
import mikera.vectorz.Vector;
import mikera.vectorz.util.VectorzException;

/**
 * Indexed sparse vector. Mutable only in the elements included in the index.
 * 
 * Index must be distinct and sorted.
 * 
 * @author Mike
 *
 */
public class SparseIndexedVector extends ASparseVector {
	private static final long serialVersionUID = 750093598603613879L;

	private final int length;
	public final Index index;
	public final double[] data;
	
	
	private SparseIndexedVector(int length, Index index) {
		this(length,index,new double[index.length()]);
	}
	
	private SparseIndexedVector(int length, Index index, double[] data) {
		this.length=length;
		this.index=index;
		this.data=data;
	}
	
	private SparseIndexedVector(int length, Index index, AVector data) {
		this.length=length;
		this.index=index;
		this.data=new double[index.length()];
		data.copyTo(this.data, 0);
	}
	
	/**
	 * Creates a SparseIndexedVector with the specified index and data values.
	 * Performs no checking - Index must be distinct and sorted.
	 */
	public static SparseIndexedVector wrap(int length, Index index, double[] data) {
		assert(index.length()==data.length);
		assert(index.isDistinctSorted());
		return new SparseIndexedVector(length, index,data);
	}
	
	public static SparseIndexedVector create(int length, Index index, double[] data) {
		if (!index.isDistinctSorted()) {
			throw new VectorzException("Index must be sorted and distinct");
		}
		if (!(index.length()==data.length)) {
			throw new VectorzException("Length of index: mismatch woth data");			
		}
		return new SparseIndexedVector(length, index,data);
	}
	
	public static SparseIndexedVector create(int length, Index index, AVector data) {
		SparseIndexedVector sv= create(length, index, new double[index.length()]);
		data.copyTo(sv.data, 0);
		return sv;
	}
	
	@Override
	public int nonSparseElementCount() {
		return data.length;
	}
	
	@Override
	public int length() {
		return length;
	}
	
	@Override
	public void multiply (AVector v) {
		if (v instanceof ArrayVector) {
			multiply((ArrayVector)v);
			return;
		}
		for (int i=0; i<data.length; i++) {
			data[i]*=v.get(index.data[i]);
		}
	}
	
	public void multiply(ArrayVector v) {
		multiply(v.getArray(),v.getArrayOffset());
	}
	
	@Override
	public void multiply(double[] array, int offset) {
		for (int i=0; i<data.length; i++) {
			data[i]*=array[offset+index.data[i]];
		}
	}
	
	@Override
	public double magnitudeSquared() {
		double result=0.0;
		for (int i=0; i<data.length; i++) {
			double d=data[i];
			result+=d*d;
		}
		return result;
	}
	
	@Override
	public boolean isZeroVector() {
		for (int i=0; i<data.length; i++) {
			if (data[i]!=0.0) return false;
		}
		return true;
	}
	
	@Override
	public double maxAbsElement() {
		double result=0.0;
		for (int i=0; i<data.length; i++) {
			double d=Math.abs(data[i]);
			if (d>result) result=d; 
		}
		return result;
	}
	
	@Override
	public void negate() {
		for (int i=0; i<data.length; i++) {
			data[i]=-data[i]; 
		}
	}
	
	@Override
	public void absolute() {
		for (int i=0; i<data.length; i++) {
			data[i]=Math.abs(data[i]); 
		}
	}

	@Override
	public double get(int i) {
		int ip=index.indexPosition(i);
		if (ip<0) return 0.0;
		return data[ip];
	}
	
	@Override
	public boolean isFullyMutable() {
		return false;
	}
	
	@Override
	public double elementSum() {
		double result=0.0;
		for (int i=0; i<data.length; i++) {
			result+=data[i];
		}
		return result;
	}
	
	
	@Override
	public double dotProduct(AVector v) {
		double result=0.0;
		for (int i=0; i<data.length; i++) {
			result+=data[i]*v.get(index.data[i]);
		}
		return result;
	}
	
	public double dotProduct(ArrayVector v) {
		double[] array=v.getArray();
		int offset=v.getArrayOffset();
		double result=0.0;
		for (int i=0; i<data.length; i++) {
			result+=data[i]*array[offset+index.data[i]];
		}
		return result;
	}
	
	@Override
	public void addMultipleToArray(double factor,int offset, double[] array, int arrayOffset, int length) {
		int start=index.seekPosition(offset);
		for (int i=start; i<data.length; i++) {
			int di=index.data[i];
			if (di>=(offset+length)) return;
			array[di+arrayOffset]+=factor*data[i];
		}
	}
	
	@Override
	public void addToArray(int offset, double[] array, int arrayOffset, int length) {
		int start=index.seekPosition(offset);
		for (int i=start; i<data.length; i++) {
			int di=index.data[i];
			if (di>=(offset+length)) return;
			array[di+arrayOffset]+=data[i];
		}
	}
	
	@Override public void copyTo(double[] array, int offset) {
		Arrays.fill(array,offset,offset+length,0.0);
		copySparseValuesTo(array,offset);
	}
	
	public void copySparseValuesTo(double[] array, int offset) {
		for (int i=0; i<data.length; i++) {
			array[offset+index.data[i]]=data[i];
		}	
	}
	
	@Override public void copyTo(AVector v, int offset) {
		if (v instanceof ArrayVector) {
			ArrayVector av=(ArrayVector)v;
			copyTo(av.getArray(),av.getArrayOffset()+offset);
		}
		v.fillRange(offset,length,0.0);
		for (int i=0; i<data.length; i++) {
			v.set(offset+index.data[i],data[i]);
		}	
	}

	@Override
	public void set(int i, double value) {
		int ip=index.indexPosition(i);
		if (ip<0) {
			throw new VectorzException("Can't set SparseIndexedVector at non-indexed position: "+i);
		}
		data[ip]=value;
	}

	@Override
	public Vector nonSparseValues() {
		return Vector.wrap(data);
	}

	@Override
	public boolean includesIndex(int i) {
		return index.indexPosition(i)>=0;
	}
	
	@Override
	public SparseIndexedVector exactClone() {
		return new SparseIndexedVector(length,index.clone(),data.clone());
	}

}
