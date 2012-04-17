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
public class Tree {
    /**************************Fields**********************************/
    private int hashval1, hashval2, noTax, taxID;
    private Vector children;
    private Tree parent;

    /*************************Constructors****************************/
    public Tree(int noTax) {
	children = new Vector();
	parent = null;
	this.noTax = noTax + 1;//Faked root with the real consensus tree as the first child and a faked node as the second child
	taxID = -1;
    }

    public Tree(BinTree bt) {
	hashval1 = bt.getHash(1);
	hashval2 = bt.getHash(2);
	noTax = bt.getNoTax();
	taxID = bt.getTaxID();
	children = new Vector();
	parent = null;
    }
    /************************Methods*********************************/
    public Tree getParent() {
	return parent;
    }

    public Vector getChildren() {
	return children;
    }

    public int getNoTax() {
	return noTax;
    }

    public int getTaxID() {
	return taxID;
    }

    public int getHashVal(int i) {
	if (i == 1) {
	    return hashval1;
	} else {
	    return hashval2;
	}
    }

    public Tree addChild (BinTree bt) {
	Tree t = new Tree(bt);
	children.add(t);
	t.setParent(this);
	return t;
    }
    public void addChild (Tree t) {
	if (!children.contains(t)) {
	    children.add(t);
	}
    }

    public void setParent(Tree t) {
	parent = t;
    }

    public void checkTree() throws DoubleCollision {
	if (children.isEmpty()) {
	    if (noTax != 1) {
		throw new DoubleCollision();
	    }
	} else {
	    Tree t;
	    Enumeration e = children.elements();
	    int sum = 0;
	    while (e.hasMoreElements()) {
		t = (Tree)e.nextElement();
		t.checkTree();
		sum += t.getNoTax();
	    }
	    if (sum != noTax) {
		throw new DoubleCollision();
	    }
	}
    }

    public void removeChild(Tree t) {
	children.remove(t);
    }

    public String toString(String[] taxa) {
	if (noTax == 1) {
	    return(taxa[taxID]);
	    //return("" + taxID);
	} else {
	    String result = "(";
	    boolean first = true;
	    Enumeration e = children.elements();
	    while (e.hasMoreElements()) {
		if (first) {
		    result += ((Tree)e.nextElement()).toString(taxa);
		    first = false;
		} else {
		    result += ("," + ((Tree)e.nextElement()).toString(taxa));
		}
	    }
	    result += ")"; //("No=" + noTax + ")");
	    return result;
	}
    }
} 
	
