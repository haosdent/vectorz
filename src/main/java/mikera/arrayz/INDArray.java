package mikera.arrayz;

import java.util.List;

import mikera.vectorz.AVector;
import mikera.vectorz.IOp;
import mikera.vectorz.Op;

/**
 * Interface for general multi-dimensional arrays of doubles
 * @author Mike
 */
public interface INDArray extends Cloneable {
	
	public int dimensionality();
	
	/**
	 * Returns the shape of the array as an array of ints.
	 * @return
	 */
	public int[] getShape();
	public int getShape(int dim);
	
	/**
	 * Returns the shape of the array as an array of longs.
	 * @return
	 */
	public long[] getLongShape();
	
	public double get();
	public double get(int x);
	public double get(int x, int y);
	public double get(int... indexes);

	public void set(double value);
	public void set(int x, double value);
	public void set(int x, int y, double value);
	public void set(int[] indexes, double value);
	public void set(INDArray a);
	public void set(Object o);
	
	public void add(INDArray a);
	public void sub(INDArray a);

	public INDArray innerProduct(INDArray a);
	public INDArray outerProduct(INDArray a);
	
	/**
	 * Creates a view of the array as a single vector in row-major order.
	 * @return
	 */
	public AVector asVector();
	
	public INDArray reshape(int... dimensions);

	public INDArray broadcast(int... dimensions);

	
	public INDArray slice(int majorSlice);
	
	public int sliceCount();
	
	public long elementCount();
	public long nonZeroCount();
	
	/**
	 * Returns true if the INDArray is mutable (at least partially)
	 * @return
	 */
	public boolean isMutable();
	
	/**
	 * Returns true if the INDArray is fully mutable in all positions
	 * i.e. every position can store any valid double value
	 * @return
	 */
	public boolean isFullyMutable();
	
	/**
	 * Returns true if the IND has additional constraints on element values
	 */
	public boolean isElementConstrained();
	
	/**
	 * Return true if this is a view
	 * @return
	 */
	public boolean isView();

	/**
	 * Returns a clone of the array, as a new array which will be fully mutable
	 * and may be of a different class to the original.
	 * @return
	 */
	public INDArray clone();

	/**
	 * Applies a unary operator to all elements of the array (in-place)
	 * @param op
	 */
	void applyOp(Op op);

	/**
	 * Applies a unary operator to all elements of the array (in-place)
	 * @param op
	 */
	void applyOp(IOp op);
	
	/**
	 * Returns true if the two arrays are exactly equal in value and shape
	 * @param a
	 * @return
	 */
	public boolean equals(INDArray a);

	/**
	 * Returns an exact deep clone of an array (i.e. of the same class as the original).
	 * @return
	 */
	public INDArray exactClone();
	
	/**
	 * Sets all elements in an array using the given double values
	 */
	public void setElements(double[] values);
	public void setElements(double[] values, int offset, int length);

	/**
	 * Gets all elements of the array, copying them into a double array
	 * @param d
	 */
	public void getElements(double[] dest, int offset);
	
	public void scale(double d);
	public void multiply(double d);

	public List<?> getSlices();
	
	/**
	 * Validates the internal data structure of the INDArray. Throws an exception on failure.
	 * 
	 * Failure indicates a serious bug and/or data corruption.
	 */
	public void validate();

	/**
	 * Copies the elements of this INDArray to the specified double array
	 * @param arr
	 */
	public void copyTo(double[] arr);
}
