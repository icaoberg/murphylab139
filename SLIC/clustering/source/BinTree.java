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
import java.io.*;
/**
 * Protocol for a binary tree
 */
public class BinTree {
    private int noTax;//Number of taxa in the tree
    private int taxID;//if leaf, the ID of the taxa
    private BinTree left, right;//two children
    private int hashval1, hashval2; //Stored two hash values

    public BinTree(String text, int hash1, int hash2, int[] weight) {
	taxID = 0;
	String s = text.trim();
	if (s.startsWith("(")) {
	    int j = 1;
	    if (s.charAt(j) == '(') {
		int t = 1;
		while (t > 0) {
		    switch (s.charAt(++j)) {
		       case '(': t++; break;
		       case ')': t--; break;
		    }
		}
		left = new BinTree(s.substring(1, j + 1), 
				   hash1, hash2, weight);
		s = s.substring(j + 1, s.length() - 1).split(",", 2)[1];
		right = new BinTree(s.substring(0, s.length()),
				    hash1, hash2, weight);
	    } else {
		String[] st = s.substring(1, s.length() - 1).split(",", 2);
		left = new BinTree(st[0], hash1, hash2, weight);
		right = new BinTree(st[1], hash1, hash2, weight);
	    }
	    noTax = left.getNoTax() + right.getNoTax();
	    hashval1 = (left.getHash(1) + right.getHash(1)) % hash1;
	    hashval2 = (left.getHash(2) + right.getHash(2)) % hash2;
	} else {
	    taxID = Integer.parseInt(s);
	    left = null;
	    right = null;
	    noTax = 1;
	    hashval1 = weight[taxID] % hash1;
	    hashval2 = weight[taxID] % hash2;
	}
    }

    public int getTaxID() {
	return taxID;
    }

    public int getNoTax() {
	return noTax;
    }

    public int getHash(int i){
	if (i == 1) {
	    return hashval1;
	}else{
	    return hashval2;
	}
    }
    public BinTree getLeft() {
	return left;
    }
    public BinTree getRight() {
	return right;
    }
    public Vector posttranverse() {
	Vector result;
	if (left == null) {
	    result = new Vector();
	} else {
	    result = left.posttranverse();
	}
	if (right != null) {
	    result.addAll(right.posttranverse());
	}
	result.add(this);
	return result;
    }
    public Vector pretranverse() {
	Vector result = new Vector();
	result.add(this);
	if (left != null) {
	    result = left.pretranverse();
	}
	if (right != null) {
	    result.addAll(right.pretranverse());
	}
	return result;
    }
}
		
