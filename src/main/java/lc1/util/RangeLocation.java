package lc1.util;

public class RangeLocation {
     int start,end;
	public RangeLocation(int start, int end) {
		this.start = start;
		this.end =end;
		// TODO Auto-generated constructor stub
	}

	
	public double overlap(RangeLocation loc){
		   return Math.min(end - loc.start , loc.end - start );
		 }
	public boolean overlaps(RangeLocation loc) {
		return overlap(loc)>0;
	}
	
	public RangeLocation intersection(RangeLocation loc) {
		int st1 = Math.max(loc.start, start);
		int end1 = Math.min(loc.end, end);
		if(st1 <= end1) return new RangeLocation(st1,end1);
		else return null;
	}
	public Location union(RangeLocation loc) {
		double overlap = overlap(loc);
		if(overlap>=0){	
			int st1 = Math.min(loc.start, start);
			int end1 = Math.max(loc.end, end);
			return new Location(st1,end1);
		}else{
			return new Location(new RangeLocation[]{this,loc});
		}
	
		
	}


	public boolean contains(int pos) {
		return pos>=start && pos < end;
	}
}
