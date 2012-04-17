/*
 * Copyright (C) 2006 Murphy Lab,Carnegie Mellon University
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation; either version 2 of the License,
 * or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 * 
 * For additional information visit http://murphylab.web.cmu.edu or
 * send email to murphy@cmu.edu
 */
import java.util.*;

public class Node{
    /***************************Fields**************************************/
    private int hashval2, noTax, taxID, count;
    private Tree tree; //used only in reconstruction
    private Node next;

    /**************************Constructors********************************/
    public Node(int hashval2, int noTax) {
	this(hashval2, noTax, -1);
    }

    public Node(int hashval2, int noTax, int taxID) {
	this.hashval2 = hashval2;
	this.noTax = noTax;
	this.taxID = taxID;
	count = 1;
	tree = null;
	next = null;
    }

    public Node(int hashval2, int noTax, int taxID, Tree t) {
	this.hashval2 = hashval2;
	this.noTax = noTax;
	this.taxID = taxID;
	count = 1;
	tree = t;
	next = null;
    }


    /*************************Methods*************************************/
    public int getHashVal() {
	return hashval2;
    }

    public int getNoTax() {
	return noTax;
    }

    public int getID() {
	return taxID;
    }

    public int getCount() {
	return count;
    }

    public Node getNext() {
	return next;
    }
    
    public void setNext(Node u) {
	next = u;
    }

    public void update(int hashval, int noTax, int ID) throws DoubleCollision{
	if ((hashval == hashval2) && (noTax == this.noTax)){
	    if (ID == taxID) {
		count++;
	    } else {
		throw new DoubleCollision();
	    }
	} else if (next != null) {
	    next.update(hashval, noTax, ID);
	} else {
	    next = new Node(hashval, noTax, ID);
	}
    }

    public void add(Tree t) {
	Node u = new Node(t.getHashVal(2), t.getNoTax(), t.getTaxID(), t);
	u.setNext(next);
	next = u;
    }
	

    public boolean check(int hashval, int noTax, double thresh){
	if (hashval == hashval2 && noTax == this.noTax) {
	    if (count > thresh) {
		return true;
	    } else {
		return false;
	    }
	} else if (next != null) {
	    return next.check(hashval, noTax, thresh);
	} else {
	    return false;
	}
    }

    public Node get(int hashval, int noTax) {
	if (hashval == hashval2 && noTax == this.noTax) {
	    return this;
	} else if (next != null) {
	    return next.get(hashval, noTax);
	} else {
	    return null;
	}
    }

    public Tree getTree(int hashval, int noTax){
	Node u = get(hashval, noTax);
	return u.tree;
    }
}
