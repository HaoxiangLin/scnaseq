package lc1.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Location {
 List<RangeLocation> locs;
 
 
 
 public Location(int start, int end){
	locs = Arrays.asList(new RangeLocation[] {new RangeLocation(start, end)});
 }
 public Location(RangeLocation[] rangeLocations) {
	this.locs = Arrays.asList(rangeLocations);
}

 
public Location intersection(Location loc) {
	List<RangeLocation> intersect = new ArrayList<RangeLocation>();
	for(int i =0; i<locs.size(); i++){
		RangeLocation loci = locs.get(i);
		for(int j=0; j<loc.locs.size(); j++){
			RangeLocation locij= loc.locs.get(j);
			if(locij.overlaps(loci)){ intersect.add(locij.intersection(loci));
			}
		}
	}
	return new Location(intersect.toArray(new RangeLocation[0]));
}
public Location union(Location loc) {
	List<RangeLocation> union= new ArrayList<RangeLocation>();
	for(int i =0; i<locs.size(); i++){
		union.add(locs.get(i));
	}
	for(int j =0; j<loc.locs.size(); j++){
		union.add(loc.locs.get(j));
	}
	Location loc1 = new Location(union.toArray(new RangeLocation[0]));
	loc1.merge();
	return loc1;
}

public void merge(){
	throw new RuntimeException ("!! need to implement");
}
public boolean overlaps(Location loc) {
	for(int i=0; i<this.locs.size(); i++){
		RangeLocation rl = locs.get(i);
		for(int j=0; j<loc.locs.size(); j++){
			if(rl.overlaps(loc.locs.get(j))) return true;
		}
	}
	return false;
}
public boolean contains(int pos) {
	for(int i=0; i<this.locs.size(); i++){
		if(locs.get(i).contains(pos)) return true;
	}
	return false;
}
}
