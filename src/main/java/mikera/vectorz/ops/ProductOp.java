package mikera.vectorz.ops;

import mikera.vectorz.Op;

public class ProductOp extends Op {
	public final Op a;
	public final Op b;
	
	private ProductOp(Op a, Op b) {
		this.a=a;
		this.b=b;
	}
	
	private static Op tryOptimisedCreate(Op a, Op b) {
		if (a instanceof ConstantOp) {
			return LinearOp.create(((ConstantOp) a).value,0.0).compose(b);
		}
		if ((a instanceof ALinearOp)&&(b instanceof ALinearOp)) {
			ALinearOp la=(ALinearOp)a;
			ALinearOp lb=(ALinearOp)b;
			double a1=la.getFactor();
			double a2=la.getConstant();
			double b1=lb.getFactor();
			double b2=lb.getConstant();		
			return QuadraticOp.create(a1*b1,a1*b2+b1*a2,b2*a2);
		}
		return null;
	}
	
	public static Op create(Op a, Op b) {
		Op t1=tryOptimisedCreate(a,b);
		if (t1!=null) return t1;
		
		Op t2=tryOptimisedCreate(b,a);
		if (t2!=null) return t2;
		
		return new ProductOp(a,b);
	}
	
	@Override
	public double apply(double x) {
		return a.apply(x)*b.apply(x);
	}

	@Override
	public double averageValue() {
		return a.averageValue()*b.averageValue();
	}
	
	@Override
	public boolean hasDerivative() {
		return (a.hasDerivative())&&(b.hasDerivative());
	}
	
	@Override
	public double derivative(double x) {
		double ay=a.apply(x);
		double by=b.apply(x);
		return a.derivative(x)*by+b.derivative(x)*ay;
	}
	
	@Override
	public Op getDerivativeOp() {
		Op da=a.getDerivativeOp();
		Op db=b.getDerivativeOp();
		return (da.product(b)).sum(db.product(a));
	}
}