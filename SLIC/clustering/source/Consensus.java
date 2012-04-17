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
public class Consensus {
    private int hash1, hash2 = 74101, noTax, noTrees;
    private int[] primes = {5, 13, 23, 31, 41, 53, 61, 71, 83, 97, 101, 113, 127, 131, 149, 151, 163, 179, 191, 211, 229, 241, 257, 271, 281, 307, 317, 347, 359, 379, 397, 419, 433, 449, 461, 479, 499, 521, 547, 571, 601, 613, 631, 641, 653, 673, 683, 701, 727, 751, 773, 797, 821, 853, 877, 911, 937, 967, 997, 1021, 1051, 1087, 1109, 1129, 1153, 1181, 1201, 1229, 1259, 1289, 1307, 1327, 1361, 1399, 1429, 1459, 1583, 1613, 1637, 1667, 1697, 1733, 1777, 1811, 1847, 1877, 1907, 1933, 1973, 2003};
    private String[] taxa; //Taxa labels
    private Vector trees; //Input trees
    private Tree consensus; //Output
    private Node[] hashTable, built;
    private double thresh;

    /**********************Constructor************************************/
    /**
     * public consensus(String filename, double perc)
     * @param String filename - the input file created by xc_bootstrap_half.m
     * @param double perc - threshold for consensus
     *
     * This constructor read in the taxa labes and the input trees and 
     * initialize the weight table and hashtables.
     */
    public Consensus(String filename, double perc) {
	BufferedReader input = null;
	// Read in the taxa information and input trees
	try {
	    input = new BufferedReader(new FileReader(filename));
	    String line = input.readLine();
	    StringTokenizer st;
	    System.out.println("Reading taxa blocks.");
	    //Read in the taxa labels
	    while (line.indexOf("dimensions ntax") == -1)
		line = input.readLine();
	    //No. of taxa
	    noTax = Integer.parseInt(line.substring(
			    line.indexOf("=") + 1, line.indexOf(";")).trim());
	    System.out.println("No. of taxa: " + noTax);
	    taxa = new String[noTax];
	    //Random weight, used in hashcode calculation
	    int[] weight = new int[noTax];
	    while (line.indexOf("taxlabels") == -1)
		line = input.readLine();
	    st = new StringTokenizer(line.substring(0, line.indexOf(";")));
	    st.nextToken();
	    int count = 0;
	    Random r = new Random();
	    while (st.hasMoreTokens()){ 
		taxa[count] = st.nextToken();
		weight[count++] = Math.abs(r.nextInt());
	    }
	    //Determine hash1
	    int low = 0;
	    int high = primes.length;
	    while (low + 1 < high) {
		//System.out.println(primes[(high + low) / 2] + " High = " +
		//		   high + " Low = " + low);
		if (primes[(high + low) / 2] < (10 * noTax)) {
		    low = (high + low) / 2;
		} else {
		    high = (high + low) / 2;
		}
	    }
	    //change from low+1 to low to avoie array index out of bound when low=length-1
	    hash1 = primes[low];
	    System.out.println("Hash1 = " + hash1);
	    hashTable = new Node[hash1];
	    built = new Node[hash1];

	    //Read in the input trees
	    System.out.print("Start reading input trees");
	    trees = new Vector();
	    while (line != null) {
		if (line.startsWith("tree B_")) {
		    trees.add(new BinTree(convertToID(line.substring(
						      line.indexOf("("), 
						      line.indexOf(";"))),
					  hash1, hash2, weight));
		    System.out.print(".");
		}
		line = input.readLine();
	    }
	    trees.trimToSize();
	    noTrees = trees.size();
	    System.out.println("\nNo. of trees read: " + noTrees);
	    thresh = noTrees * perc;
	} catch (Exception e) {
		System.out.println("Exception: " + e);
		e.printStackTrace();
		System.exit(1);
	} finally {
	    try {
		input.close();
	    } catch (IOException ie) {
		System.out.println("Exception: " + ie);
		ie.printStackTrace();
		System.exit(1);
	    }
	}
    }
	
    /**************************Methods*************************************/
    /**
     * private String convertToID(String originalTree)
     * @param String originalTree - The input tree directly read from file
     * @return A String where the taxa labels are replaced by their
     * corresponding index in the array taxa
     */
    private String convertToID(String originalTree) {
	String result = new String(originalTree);
	for (int i = 0; i < noTax; i++) {
	    int indx = result.indexOf(taxa[i] + ")");
	    if (indx == -1) {
		indx = result.indexOf(taxa[i] + ",");
	    }
	    try{
		result = result.substring(0, indx) + i + 
		    result.substring(indx + taxa[i].length(), result.length());
	    } catch (Exception e) {
		System.out.println("Error in replacing\n" + taxa[i] +
				   "\n in " + result);
		System.exit(-1);
	    }
	}
	return result;
    }
		    
    public void process() throws DoubleCollision {
	System.out.print("Start processing the input trees");
	Enumeration e = trees.elements();
	while (e.hasMoreElements()) {
	    BinTree bt = (BinTree)e.nextElement();
	    Enumeration eb = bt.posttranverse().elements();
	    while (eb.hasMoreElements()) {
		BinTree tmp = (BinTree)eb.nextElement();
		try {
		    if (hashTable[tmp.getHash(1)] != null) {
			hashTable[tmp.getHash(1)].update(tmp.getHash(2),
							 tmp.getNoTax(),
							 tmp.getTaxID());
		    } else {
			hashTable[tmp.getHash(1)]= new Node(tmp.getHash(2),
							    tmp.getNoTax(),
							    tmp.getTaxID());
		    }
		} catch (ArrayIndexOutOfBoundsException ae) {
		    System.out.println("\nException in process: Hashcode(" +
				       tmp.getHash(1) + "), Hash1 = " + hash1);
		    System.exit(-1);
		}
	    }
	    System.out.print(".");
	}
	System.out.println("\nFinishing processing input trees.");
    }

    private void constructFromBin(BinTree bt, Tree p) {
	Tree t;
	if (hashTable[bt.getHash(1)].check(bt.getHash(2), bt.getNoTax(), thresh)) {
	    if ((built[bt.getHash(1)] == null) ||
		(built[bt.getHash(1)].get(bt.getHash(2), bt.getNoTax()) == null)) {//Not exists
		t = p.addChild(bt);
		if (built[bt.getHash(1)] == null) {
		    built[bt.getHash(1)] = new Node(bt.getHash(2), bt.getNoTax(), bt.getTaxID(), t);
		} else {
		    built[bt.getHash(1)].add(t);
		}
	    } else { // check for parents
		t = built[bt.getHash(1)].getTree(bt.getHash(2), bt.getNoTax());
		if (t.getParent().getNoTax() > p.getNoTax()) {
		    t.getParent().removeChild(t);
		    p.addChild(t);
		    t.setParent(p);
		}
	    }
	} else {
	    t = p;
	}
	if (bt.getLeft() != null) {
	    constructFromBin(bt.getLeft(), t);
	}
	if (bt.getRight() != null) {
	    constructFromBin(bt.getRight(), t);
	}
    }

	
    public void construct() throws DoubleCollision{
	System.out.print("Start constructing the consensus tree");
	Tree result = new Tree(noTax);
	Enumeration e = trees.elements();
	while (e.hasMoreElements()) {
	    constructFromBin((BinTree)e.nextElement(), result);
	    System.out.print(".");
	}
	System.out.println("\nFinishing constructing the consensus tree.");
      	consensus =  (Tree)result.getChildren().get(0);
	//System.out.println(toString());
	consensus.checkTree();
    }

    public String toString() {
	String result = "#NEXUS\n\nbegin taxa;\ndimensions ntax = " + 
	                 consensus.getNoTax() +  ";\ntaxlabels";
	for (int i = 0; i < taxa.length; i++)
	    result += (" " + taxa[i]);
	result += ";\nend;\n\nbegin trees;\ntree *MAJORITY_RULE-COMPONENT = [&R]\n" + consensus.toString(taxa) + ";\nend;";
	return result;
    }
}
