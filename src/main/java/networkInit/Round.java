package networkInit;

//import java.io.Serializable;

public class Round {

	/**
	 * 
	 */
	private float x, y;
	private float radius;

	public Round(float x, float y, float radius) {
		//super();
		this.x = x;
		this.y = y;
		this.radius = radius;
	}

	public float getX() {
		return x;
	}

	public void setX(float x) {
		this.x = x;
	}

	public float getY() {
		return y;
	}

	public void setY(float y) {
		this.y = y;
	}

	public float getRadius() {
		return radius;
	}

	public void setRadius(float radius) {
		this.radius = radius;
	}

	@Override
	public String toString() {
		return "Round [x=" + x + ", y=" + y + ", radius=" + radius + "]";
	}

	public float getDis(Round round) {
		return (float) Math.sqrt((this.x - round.x) * (this.x - round.x) + (this.y - round.y) * (this.y - round.y));
	}

	public float getArea() {
		return (float) (Math.PI * this.radius * this.radius);
	}

}
