package mikera.vectorz;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import mikera.arrayz.INDArray;
import mikera.arrayz.TestArrays;
import mikera.indexz.Index;
import mikera.indexz.Indexz;
import mikera.matrixx.AMatrix;
import mikera.matrixx.Matrixx;
import mikera.matrixx.impl.MatrixAsVector;
import mikera.util.Rand;
import mikera.vectorz.impl.ArraySubVector;
import mikera.vectorz.impl.AxisVector;
import mikera.vectorz.impl.ImmutableVector;
import mikera.vectorz.impl.IndexVector;
import mikera.vectorz.impl.RangeVector;
import mikera.vectorz.impl.RepeatedElementVector;
import mikera.vectorz.impl.IndexedArrayVector;
import mikera.vectorz.impl.IndexedSubVector;
import mikera.vectorz.impl.JoinedArrayVector;
import mikera.vectorz.impl.SingleElementVector;
import mikera.vectorz.impl.SparseHashedVector;
import mikera.vectorz.impl.SparseImmutableVector;
import mikera.vectorz.impl.SparseIndexedVector;
import mikera.vectorz.impl.StridedVector;
import mikera.vectorz.impl.Vector0;
import mikera.vectorz.impl.WrappedSubVector;
import mikera.vectorz.ops.Constant;
import mikera.vectorz.util.VectorBuilder;

import org.junit.Test;


public class TestVectors {
	@Test public void testDistances() {
		Vector3 v=Vector3.of(1,2,3);
		Vector3 v2=Vector3.of(2,4,6);
		assertEquals(v.magnitude(),v.distance(v2),0.000001);
		assertEquals(6,v.distanceL1(v2),0.000001);
		assertEquals(3,v.distanceLinf(v2),0.000001);
	}
	
	public void testDistances(AVector v) {
		AVector zero=Vectorz.newVector(v.length());
		assertEquals(v.maxAbsElement(),v.distanceLinf(zero),0.0);
		
		if (!v.isFullyMutable()) return;
		AVector z=v.exactClone();
		z.fill(0.0);
		assertEquals(v.maxAbsElement(),v.distanceLinf(z),0.0);
	}
	
	public void testMagnitudes(AVector v) {
		double d=v.magnitude();
		double d2=v.magnitudeSquared();
		assertEquals(d*d,d2,0.00001);
		
		assertTrue(d<=(v.maxAbsElement()*v.length()));
	}
	
	public void testSquare(AVector v) {
		if (!v.isFullyMutable()) return;
		v=v.exactClone();
		AVector vc=v.clone();
		v.square();
		vc.square();
		assertEquals(vc,v);
	}
	
	@Test public void testCross() {
		Vector3 v=Vector3.of(1,2,3);
		v.crossProduct(Vector.of(1,1,1));
		assertEquals(Vector3.of(-1,2,-1),v);
	}
	
	@Test public void testToString() {
		Vector3 v=Vector3.of(1,2,3);
		assertEquals("[1.0,2.0,3.0]",v.toString());
	}
	
	@Test public void testClamping() {
		Vector3 v=Vector3.of(1,2,3);
		v.clamp(1.5, 2.5);
		assertEquals(Vector3.of(1.5,2,2.5),v);
	}
	
	@Test public void testClampMin() {
		Vector3 v=Vector3.of(1,2,3);
		v.clampMin(1.5);
		assertEquals(Vector3.of(1.5,2,3),v);
	}
	
	@Test public void testClampMax() {
		Vector3 v=Vector3.of(1,2,3);
		v.clampMax(2.5);
		assertEquals(Vector3.of(1,2,2.5),v);
	}
	
	@Test public void testElementSum() {
		Vector3 v=Vector3.of(1,2,3);
		assertEquals(6.0,v.elementSum(),0.0);
	}
	
	
	@Test public void testCreateFromIterable() {
		ArrayList<Object> al=new ArrayList<Object>();
		al.add(1);
		al.add(2L);
		al.add(3.0);
		AVector v=Vectorz.create((Iterable<Object>)al);
		assertEquals(Vector.of(1,2,3),v);
	}
	
	@Test public void testSubVectors() {
		double[] data=new double[100];
		for (int i=0; i<100; i++) data[i]=i;
		
		ArraySubVector v=ArraySubVector.wrap(data);
		assertEquals(10,v.get(10),0.0);
		assertTrue(v.isView());
		
		AVector v2=v.subVector(5, 90);
		assertEquals(90,v2.length());
		assertEquals(15,v2.get(10),0.0);
		assertTrue(v2.isView());
		
		AVector v3=v2.subVector(5,80);
		assertEquals(20,v3.get(10),0.0);
		assertTrue(v3.isView());
		
		v3.set(10, -99);
		assertEquals(-99,v.get(20),0.0);
	}
	
	@Test public void testWrap() {
		double[] data=new double[]{1,2,3};
		
		Vector v=Vector.wrap(data);
		data[0]=13;
		assertEquals(13,v.get(0),0.0);
	}
	
	@Test public void testSubVectorCopy() {
		double[] data=new double[100];
		for (int i=0; i<100; i++) data[i]=i;
		ArraySubVector v=ArraySubVector.wrap(data);	
		
		assertEquals(Arrays.hashCode(data),v.hashCode());
		
		ArraySubVector v2=v.exactClone();
		assertEquals(v,v2);
		assertTrue(v2.isView());
		
		v.set(10,Math.PI);
		assertTrue(!v.equals(v2));
	}
	
	@Test public void testString() {
		assertEquals("[0.0,0.0]",new Vector2().toString());
		assertEquals("[1.0,2.0]",Vectorz.create(1,2).toString());
	}
	
	private void testParse(AVector v) {
		String str = v.toString();
		AVector v2=Vectorz.parse(str);
		if (v.length()>0) {
			assertEquals(v.get(0),v2.get(0),0.0);			
		}
		assertEquals(v,v2);
	}
	
	private void testHashCode(AVector v) {
		assertEquals(v.hashCode(),v.toList().hashCode());
	}
	
	private void testAddToArray(AVector v) {
		int len=v.length();
		double[] ds=new double[len+10];
		if (len>=3) {
			v.addToArray(2, ds, 5, len-2);
			assertEquals(v.get(2),ds[5],0.000);
			assertEquals(0.0,ds[5+len],0.0);
		}
		
		v.addToArray(0, ds, 0, len);
		if (len>0) {
			assertEquals(ds[0],v.get(0),0.0);
		}
		
		double[] ds2=new double[len];
		v.addToArray(ds2, 0);
		assertEquals(v,Vector.wrap(ds2));

	}
	
	private void testAddProduct(AVector v) {
		int len=v.length();
		if (!v.isFullyMutable()) return;
		
		v=v.exactClone();
		AVector v2=v.exactClone();
		AVector vc=new Vector(v);
		
		AVector p1=Vectorz.createUniformRandomVector(len);
		AVector p2=Vectorz.createUniformRandomVector(len);
		
		v.addProduct(p1, p2,3.0);
		vc.addProduct(p1, p2, 3.0);
		for (int i=0; i<3; i++) v2.addProduct(p1, p2);
		
		assertTrue(v.epsilonEquals(vc));
		assertTrue(v.epsilonEquals(v2));
	}
	
	private void testInnerProducts(AVector v) {
		int len=v.length();
		AVector c=Vectorz.createUniformRandomVector(v.length());
		assertEquals(v.dotProduct(c),v.innerProduct((INDArray)c).get(),0.00001);
		
		if (len>20) return;
		
		AMatrix m=Matrixx.createRandomMatrix(len, len);
		assertEquals(v.innerProduct(m),m.getTranspose().transform(v));
	}
	
	private void testAddMultipleToArray(AVector v) {
		int len=v.length();
		double[] ds=new double[len+10];
		if (len>=3) {
			v.addMultipleToArray(3.0,2, ds, 5, len-2);
			assertEquals(v.get(2)*3.0,ds[5],0.000);
			assertEquals(0.0,ds[5+len],0.0);
		}
		
		v.addMultipleToArray(2.0,0, ds, 0, len);
		if (len>0) {
			assertEquals(v.get(0)*2.0,ds[0],0.0);
		}
	}
	
	private void testAddMultiple(AVector v) {
		int len=v.length();
		Vector tv=Vector.createLength(len+10);
		tv.addMultiple(5,v,10.0);
		assertEquals(0.0,tv.get(4),0.0);
		assertEquals(10.0*v.get(0),tv.get(5),0.0);
		assertEquals(10.0*v.get(len-1),tv.get(5+len-1),0.0);
		assertEquals(0.0,tv.get(5+len),0.0);
	}
	
	private void testElementSum(AVector v) {
		double res=0.0;
		for (int i=0; i<v.length(); i++) {
			res+=v.get(i);
		}
		assertEquals(res,v.elementSum(),0.00001);
	}
	
	private void testElementProduct(AVector v) {
		double res=1.0;
		for (int i=0; i<v.length(); i++) {
			res*=v.get(i);
		}
		assertEquals(res,v.elementProduct(),0.00001);
	}
	
	
	private void testMinMax(AVector v) {
		if (v.length()==0) return;
		assertEquals(v.maxAbsElement(),Math.abs(v.get(v.maxAbsElementIndex())),0.0);
		assertEquals(v.maxElement(),v.get(v.maxElementIndex()),0.0);
		assertEquals(v.minElement(),v.get(v.minElementIndex()),0.0);
	}
	
	private void testAddMultipleIndexed(AVector v) {
		v=v.exactClone();
		if (v.isFullyMutable()) {
			Vectorz.fillGaussian(v);
		}
		int len=v.length();
		int tlen=len+10;
		
		Index ix=Indexz.createRandomChoice(len, tlen);
		Vector tv=Vector.createLength(tlen);
		
		tv.addMultiple(v, ix, 2.0);
		assertEquals(v.elementSum()*2.0,tv.elementSum(),0.0001);
	}
	
	private void testAddMultipleIndexed2(AVector v) {
		if (v.length()<1) return;
		if (!v.isFullyMutable()) return;
		v.fill(1.0);
		v=v.exactClone();
		int len=v.length();
		
		Vector tv=Vector.createLength(len);
		Vectorz.fillGaussian(tv);
		Index ix=Indexz.createLength(len);
		for (int i=0; i<len; i++) {
			ix.set(i, Rand.r(len));
		}
		
		v.addMultiple(tv, ix, 2.0);
		assertEquals(len+tv.elementSum()*2.0,v.elementSum(),0.0001);
	}
	

	private void testSubvectors(AVector v) {
		int len=v.length();
		assertEquals(v,v.subVector(0, len));
		assertEquals(0,v.subVector(0, 0).length());
		assertEquals(0,v.subVector(len, 0).length());
		int m=len/2;
		
		AVector v2=v.subVector(0, m).join(v.subVector(m, len-m));
		assertEquals(v,v2);
		
		AVector zv=v.subVector(len/2, 0);
		assertTrue(zv==Vector0.INSTANCE);
		
		// whole subvector should return same vector
		if (len>0) {
			AVector fv=v.subVector(0, len);
			assertTrue(fv==v);
		}
	}
	
	private void testClone(AVector v) {
		AVector cv=v.clone();
		int len=cv.length();
		assertEquals(v.length(), len);
		assertFalse(cv.isView());
		assertTrue((cv.length()==0)||cv.isMutable());		
		assertTrue(cv.isFullyMutable());
		assertEquals(v,cv);
		
		for (int i=0; i<len; i++) {
			double x=cv.get(i);
			// clone should have equal values
			assertEquals(v.get(i),x,0.0000001);
			cv.set(i, x+1);
			
			// changing clone should not change original
			assertNotSame(v.get(i),cv.get(i));
		}
	}
	
	private void testExactClone(AVector v) {
		AVector cv=v.exactClone();
		AVector cv2=v.clone();
		assertEquals(v,cv);
		assertEquals(v.getClass(),cv.getClass());
		
		// test that trashing the exact clone doesn't affect the original vector
		if (cv.isFullyMutable()) {
			cv.fill(Double.NaN);
			assertEquals(cv2,v);
		}
	}
	
	private void testSet(AVector v) {
		v= (v.isFullyMutable()) ? v.exactClone() : v.clone();
		
		int len=v.length();
		
		for (int i=0; i<len; i++) {
			v.set(i,i+0.5);
			assertEquals(i+0.5,v.get(i),0.0);
		}
	}
	
	private void testImmutable(AVector v) {
		AVector iv=v.immutable();
		assertFalse(iv.isMutable());
		
		try {
			iv.set(0,1.0);
			fail();
		} catch (Throwable t) {
			// OK
		}
	}
	
	
	private void testNonZero(AVector v) {
		int len=v.length();
		int n=(int)v.nonZeroCount();
		int[] nzi = v.nonZeroIndices();
		Index ind=Index.of(nzi);
		assertEquals(n,nzi.length);
		
		for (int i=0; i<nzi.length; i++) {
			assertTrue(v.unsafeGet(nzi[i])!=0.0);
		}
		
		for (int i=0; i<len; i++) {
			if (v.unsafeGet(i)==0.0) {
				assertFalse(ind.containsSorted(i));
			} else {
				assertTrue(ind.containsSorted(i));				
			}
		}
	}
	
	private void testJoining(AVector v) {
		AVector vv=v.join(v);
		assertEquals(vv,v.join(v,0));
		
		int n=v.length();
		AVector v2=vv.subVector(n, n);
		assertEquals(v,v2);
	}
	
	private void testSetElements(AVector v) {
		if (!v.isFullyMutable()) return;
		
		v=v.exactClone();
		Vectorz.fillGaussian(v);
		AVector vc=v.clone();
		
		double[] data=v.toDoubleArray();
		v.fill(0.0);
		v.setElements(data);
		
		assertEquals(vc,v);
	}
	
	private void testAdd(AVector v) {
		v=v.exactClone();
		int len=v.length();
		int split=len/2;
		int slen=len-split;
		
		AVector t=Vectorz.newVector(slen);
		t.add(v,split);
		assertEquals(t,v.subVector(split, slen));
		
		t.multiply(0.5);
		t.set(v,split);
		assertEquals(t,v.subVector(split, slen));
		t.addProduct(t, Vectorz.createZeroVector(slen));
		assertEquals(t,v.subVector(split, slen));
	}
	
	
	private void testAddEquivalents(AVector v) {
		v=v.exactClone();
		int len=v.length();
		
		AVector c=Vectorz.newVector(len);
		c.add(v);
		assertEquals(v,c);
		
		c.fill(0.0);
		c.add(0,v);
		assertEquals(v,c);
		
		c.fill(0.0);
		c.add(v,0);
		assertEquals(v,c);
		
		c.fill(0.0);
		c.addMultiple(0,v,1.0);
		assertEquals(v,c);
		
		c.fill(0.0);
		c.addMultiple(v,1.0);
		assertEquals(v,c);
				
		if (!v.isFullyMutable()) return;
	}

	
	private void testAddFromPosition(AVector v) {
		if (!v.isFullyMutable()) return;
		AVector tv=v.exactClone();
		int len=tv.length();
		
		AVector sv=Vectorz.createRange(len+10);

		tv.add(sv,5);
		assertEquals(v.get(0)+5.0,tv.get(0),0.0001);
		assertEquals(v.get(len-1)+5.0+len-1,tv.get(len-1),0.0001);
	}
	
	private void testAddToPosition(AVector v) {
		int len=v.length();
		
		AVector tv=Vectorz.createRange(len+10);

		tv.add(5,v,0,len);
		
		assertEquals(5.0+v.get(0),tv.get(5),0.0001);
		assertEquals(5.0+len,tv.get(5+len),0.0001);
	}
	
	private void testAddAt(AVector v) {
		if (!v.isFullyMutable()) return;
		int len=v.length();
		v=v.exactClone();

		for (int i=0; i<len; i++) {
			double old=v.get(i);
			v.addAt(i, i+1);
			assertEquals(old+i+1,v.get(i),0.0000001);
		}
	}
	
	private void testOutOfBoundsSet(AVector v) {
		if (!v.isFullyMutable()) return;
		try {
			v.set(-1, 0.0);
			fail("Should be out of bounds!");
		} catch (Throwable x) {
			// OK!
		}
		
		try {
			if (v instanceof GrowableVector) return;
			v.set(v.length(), 0.0);
			fail("Should be out of bounds!");
		} catch (Throwable x) {
			// OK!
		}
	}
	
	
	private void testOutOfBoundsGet(AVector v) {
		try {
			v.get(-1);
			fail("Should be out of bounds!");
		} catch (Throwable x) {
			// OK!
		}
		
		try {
			if (v instanceof GrowableVector) return;
			v.get(v.length());
			fail("Should be out of bounds!");
		} catch (Throwable x) {
			// OK!
		}
	}
	public void testSubVectorMutability(AVector v) {
		if (!v.isFullyMutable()) return;
		v=v.exactClone();
		
		int vlen=v.length();
		
		int start=Math.min(vlen/2,vlen-1);
		
		AVector s1 = v.subVector(start, vlen-start);
		AVector s2 = v.subVector(start, vlen-start);
		
		int len=s1.length();
		for (int i=0; i<len; i++) {
			double x=s2.get(i);
			// clone should have equal values
			assertEquals(s1.get(i),x,0.0000001);
			s1.set(i, x+1);
			
			// change should be reflected in both subvectors
			assertEquals(s1.get(i),s2.get(i),0.0000001);
		}
	}
	
	public void testVectorMutability(AVector v) {
		v=v.exactClone();
		if (v.isFullyMutable()) {
			assertTrue(v.isMutable());
			for (int i=0; i<v.length(); i++) {
				v.set(i,i);
				assertEquals(i,v.get(i),0.0);
			}
		}
		
		if (v.isMutable()) {
			// v.set(0,0.01); not always, may not be fully mutable
		}
	}
	
	
	private void testZero(AVector v) {
		if (!v.isFullyMutable()) return;
		v=v.exactClone();
		v.multiply(0.0);
		assertTrue(v.isZero());
	}
	
	private void testSparseOps(AVector v) {
		assertTrue(v.nonZeroCount()<=v.nonSparseIndexes().length());
		
		AVector sc=v.sparseClone();
		assertTrue("clone = "+sc.getClass(), sc.isFullyMutable());
		assertEquals(sc,v);
	}

	private void testNormalise(AVector v) {
		if (!v.isFullyMutable()) return;
		v=v.exactClone();
		
		v.set(0,v.get(0)+Math.random());
		AVector v2=v.toNormal();

		double d=v.magnitude();
		double nresult=v.normalise();
		
		assertTrue(v2.epsilonEquals(v)); // compared normalised versions
		assertEquals(d,nresult,0.0000001);
		assertTrue(v.isUnitLengthVector());
	}
	
	private void testFilling(AVector v) {
		if (!v.isFullyMutable()) return;
		v=v.exactClone();
		v.fill(1.23);
		assertEquals(1.23,Vectorz.minValue(v),0.0);
		assertEquals(1.23,Vectorz.maxValue(v),0.0);
		
		v.fillRange(0, v.length(), 1.24);
		assertEquals(1.24,Vectorz.minValue(v),0.0);
		assertEquals(1.24,Vectorz.maxValue(v),0.0);
		assertEquals(1.24,Vectorz.averageValue(v),0.0001);
	}
	
	private void testEquality(AVector v) {
		assertEquals(v,v.clone());
		assertNotEquals(v,v.join(Vector.of(1)));
		assertEquals(v,Vectorz.createSparse(v));
	}
	
	private void testCopyTo(AVector v) {
		int len=v.length();
		Vector tv=Vector.createLength(len+2);
		tv.fill(Double.NaN);
		v.copyTo(tv, 1);
		assertTrue(Double.isNaN(tv.get(0)));
		assertTrue(Double.isNaN(tv.get(len+1)));
		assertFalse(Double.isNaN(tv.get(1)));
		assertFalse(Double.isNaN(tv.get(len)));
	}

	private void testUnsafeGet(AVector v) {
		int len=v.length();
		for (int i=0; i<len; i++) {
			assertEquals(v.get(i),v.unsafeGet(i),0.0);
		}
	}
	
	private void testOutOfBounds(AVector v) {
		try {
			v.get(-1);
			fail();
		} catch (IndexOutOfBoundsException a) {/* OK */}
		
		try {
			v.get(v.length());
			fail();
		} catch (IndexOutOfBoundsException a) {/* OK */}
		
		if (v.isFullyMutable()) {
			v=v.exactClone();
			try {
				v.set(-1,1);
				fail();
			} catch (IndexOutOfBoundsException a) {/* OK */}
			
			if (!(v instanceof GrowableVector)) try {
				v.set(v.length(),1);
				fail();
			} catch (IndexOutOfBoundsException a) {/* OK */}
		}
	}
	
	private void testUnsafeSet(AVector v) {
		if (!v.isFullyMutable()) return;
		AVector a=v.exactClone();
		AVector b=v.exactClone();
		int len=v.length();
		for (int i=0; i<len; i++) {
			v.set(i, a.get(i)+1.0);
			v.unsafeSet(i, a.get(i)+1.0);
		}
		assertEquals(a,b);
	}

	private void testAsList(AVector v) {
		List<Double> al=v.asElementList();
		List<Double> tl=v.toList();
		assertEquals(al,tl);
		
		int len=v.length();
		assertEquals(len,al.size());
		assertEquals(len,tl.size());
	}
	
	private void testMultiply(AVector v) {
		if (!v.isFullyMutable()) return;
		
		int len=v.length();
		v=v.exactClone();
		
		AVector m=Vectorz.newVector(len);
		m.fill(2);
		AVector v2=v.clone();
		v2.multiply(m);
		v.multiply(2);
		assertTrue(v.epsilonEquals(v2,0.00001));
	}
	

	
	private void testDivide(AVector v) {
		if (!v.isFullyMutable()) return;

		v=v.exactClone();
		
		Vectorz.fillGaussian(v);
		
		AVector m=v.clone();
		m.multiply(v);
		m.divide(v);
		assertTrue(v.epsilonEquals(m,0.00001));
	}
	
	private void testIterator(AVector v) {
		int count=0;
		double total=0.0;
		
		for (double d: v) {
			count++;
			total+=d;
		}
		
		assertEquals(v.length(),count);
		assertEquals(Vectorz.totalValue(v),total,0.00000001);
	}
	
	private void testApplyOp(AVector v) {
		if (!v.isFullyMutable()) return;
		AVector c=v.exactClone();
		AVector d=v.exactClone();
		
		c.fill(5.0);
		d.applyOp(Constant.create(5.0));
		assertTrue(c.equals(d));
	}
	
	private void doNonDegenerateTests(AVector v) {
		if (v.length()==0) return;
		testSubVectorMutability(v);
		testSetElements(v);
		testAddMultiple(v);
		testAddMultipleIndexed(v);
		testAddMultipleIndexed2(v);
		testAddFromPosition(v);
		testAddToPosition(v);
		testVectorMutability(v);
		testCopyTo(v);
		testNormalise(v);
		testFilling(v);
	}
	
	private void testSlicing(AVector v) {
		for (int i=0; i<v.length(); i++) {
			AScalar ss=v.slice(0, i);
			if (v.isFullyMutable()) {
				assertTrue(ss.isView());
				assertTrue(ss.isFullyMutable());
			}
		}	
	}
	
	public void doGenericTests(AVector v) {
		testClone(v);
		testExactClone(v);
		testAdd(v);
		testAddEquivalents(v);
		testAddToArray(v);
		testElementSum(v);
		testElementProduct(v);
		testAddAt(v);
		testUnsafeGet(v);
		testUnsafeSet(v);
		testOutOfBounds(v);
		testSlicing(v);
		testMinMax(v);
		testJoining(v);
		testNonZero(v);
		testAddProduct(v);
		testAddMultipleToArray(v);
		testApplyOp(v);
		testEquality(v);
		testInnerProducts(v);
		testMultiply(v);
		testDivide(v);
		testSet(v);
		testSparseOps(v);
		testSquare(v);
		testSubvectors(v);
		testParse(v);
		testDistances(v);
		testMagnitudes(v);
		testZero(v);
		testHashCode(v);
		testAsList(v);
		testIterator(v);
		testOutOfBoundsSet(v);
		testOutOfBoundsGet(v);
		testImmutable(v);
		
		doNonDegenerateTests(v);
		
		new TestArrays().testArray(v);
	}

	@Test public void genericTests() {
		doGenericTests(Vector0.of());
		
		doGenericTests(new Vector1(1.0));
		doGenericTests(new Vector2(1.0,2.0));
		doGenericTests(new Vector3(1.0,2.0,3.0));
		doGenericTests(new Vector4(1.0,2.0,3.0,4.0));
		doGenericTests(new Vector4(1.0,2.0,3.0,4.0).subVector(1, 2));
		
		// bit vectors
		doGenericTests(BitVector.of());
		doGenericTests(BitVector.of(0));
		doGenericTests(BitVector.of(0,1,0));
		doGenericTests(BitVector.of(0,1,0).subVector(1, 1));
		
		// zero-length Vectors
		doGenericTests(Vector.of());
		doGenericTests(Vector.EMPTY);
		doGenericTests(new GrowableVector(Vector.of()));
		doGenericTests(Vector.wrap(new double[0]));
		doGenericTests(new Vector3(1.0,2.0,3.0).subVector(2, 0));
		
		for (int j=0; j<10; j++) {
			double[] data=new double[j];
			for (int i=0; i<j; i++) data[i]=i;
			doGenericTests(Vectorz.create(data));
		}
		
		double[] data=new double[100];
		int[] indexes=new int[100];
		for (int i=0; i<100; i++) {
			data[i]=i;
			indexes[i]=i;
		}

		doGenericTests(ArraySubVector.wrap(data));
		doGenericTests(IndexedArrayVector.wrap(data,indexes));
		doGenericTests(IndexedSubVector.wrap(Vector.of(data),indexes));
		
		doGenericTests(Vector.create(data).subVector(25, 50));
		doGenericTests(ArraySubVector.wrap(data).subVector(25, 50));
		
		AVector v3 = new Vector3(1.0,2.0,3.0);
		doGenericTests(v3.subVector(1, 2));
		doGenericTests(new WrappedSubVector(v3,1,2));

		AVector joined = Vectorz.join(v3, Vectorz.create(data));
		doGenericTests(joined);
		
		AVector v4 = Vectorz.create(1.0,2.0,3.0,4.0);
		doGenericTests(v4);
		doGenericTests(v4.subVector(1, 2));
		
		AVector g0=new GrowableVector();
		doGenericTests(g0);
		
		AVector g4=new GrowableVector(v4);
		doGenericTests(g4);
		
		AVector j5=Vectorz.join(g4,joined,v3,v4,g0,g0,joined);
		doGenericTests(j5);
		
		AMatrix m1=Matrixx.createRandomSquareMatrix(5);
		doGenericTests(m1.asVector());
		doGenericTests(m1.getRow(4));
		doGenericTests(m1.getColumn(1));
		doGenericTests(m1.getLeadingDiagonal());
		
		AMatrix m2=Matrixx.createRandomSquareMatrix(3);
		doGenericTests(m2.asVector());
		doGenericTests(m2.getRow(1));
		doGenericTests(m2.getColumn(1));
		doGenericTests(m2.getLeadingDiagonal());
		doGenericTests(new MatrixAsVector(m2));

		AMatrix m3=Matrixx.createRandomMatrix(4,5);
		doGenericTests(m3.asVector());
		doGenericTests(m3.getRow(2));
		doGenericTests(m3.getColumn(2));
		doGenericTests(m3.subMatrix(1, 1, 2, 3).asVector());
		doGenericTests(new MatrixAsVector(m3));
		
		doGenericTests(AxisVector.create(1,3));
		doGenericTests(AxisVector.create(0,1));
		doGenericTests(AxisVector.create(5,10));
		
		doGenericTests(new SingleElementVector(1,3));
		doGenericTests(new SingleElementVector(0,1));

		doGenericTests(RepeatedElementVector.create(1,1.0));
		doGenericTests(RepeatedElementVector.create(4,0.0));
		doGenericTests(RepeatedElementVector.create(10,1.0));
		doGenericTests(RepeatedElementVector.create(10,1.0).subVector(2, 5));
		
		doGenericTests(IndexVector.of(1,2,3));
		doGenericTests(IndexVector.of(1));
		doGenericTests(IndexVector.of());
		
		doGenericTests(SparseIndexedVector.create(10,Index.of(1,3,6),Vector.of(1.0,2.0,3.0)));
		doGenericTests(SparseIndexedVector.create(10,Index.of(),Vector.of()));
		doGenericTests(SparseIndexedVector.create(Vector.of(1,2,3,4,5))); // fully dense!
		doGenericTests(SparseIndexedVector.create(Vector.of(-1,-2,-3))); // fully dense!

		doGenericTests(SparseImmutableVector.create(10,Index.of(1,3,6),Vector.of(1.0,2.0,3.0)));
		doGenericTests(SparseImmutableVector.create(Vector.of(1,2,3,4,5))); // fully dense!
		doGenericTests(SparseImmutableVector.create(Vector.of(-1,-2,-3))); // fully dense!

		doGenericTests(Vector3.of(1,2,3).join(SparseIndexedVector.create(5,Index.of(1,3),Vector.of(1.0,2.0))));
		
		doGenericTests(SparseHashedVector.create(Vector.of(0,1,-1.5,0,2)));
		doGenericTests(SparseHashedVector.create(Vector.of(1,2,3,4,5))); // fully dense!
		doGenericTests(SparseHashedVector.create(Vector.of(-1,-2,-3))); // fully dense!
		
		doGenericTests(new Scalar(1.0).asVector());
		doGenericTests(Vector.of(1,2,3).slice(1).asVector());
		
		AVector jav1=JoinedArrayVector.create(v4);
		AVector jav2=JoinedArrayVector.create(j5);
		doGenericTests(jav1);
		doGenericTests(jav2);
		doGenericTests(jav2.join(jav1));
		doGenericTests(jav2.join(jav1).subVector(2, 5));
		doGenericTests(Vector3.of(1,2,3).join(JoinedArrayVector.create(g4)));
		
		doGenericTests(StridedVector.wrap(new double[]{1,2,3}, 2, 3, -1));
		doGenericTests(StridedVector.wrap(new double[]{1,2}, 1, 1, 100));
		
		doGenericTests(ImmutableVector.create(Vector.of(1,2,3)));
		doGenericTests(ImmutableVector.create(Vector.of()));
		
		doGenericTests(RangeVector.create(-10,3));
		doGenericTests(RangeVector.create(0,7));
		
		// VectorBuilder as a Vector
		VectorBuilder vbl=new VectorBuilder();
		doGenericTests(vbl);
		vbl.append(1.0);
		doGenericTests(vbl);
		vbl.append(2.0,3.0);
		doGenericTests(vbl);
		vbl.append(4,5,6);
		doGenericTests(vbl);
		
	}
}
