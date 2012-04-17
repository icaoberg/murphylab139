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
import java.io.*;
public class MajConsensus{
    public static void main(String[] argv) {
	boolean done = false;
	PrintWriter pw = null;
	String outputfile = null;
	switch (argv.length) {
	    case 1:
		outputfile = argv[0] + "_MRC";
		break;
	    case 2:
		outputfile = argv[1];
		break;
	    default:
		System.out.println("Error in using MajConsensus.\nUsage: java MajConsensus inputfile [outputfile]");
		System.exit(-1);
	}
	try {
	    pw = new PrintWriter(new FileWriter(outputfile));
	    if (pw != null) {
		while (!done) {
		    try {
			Consensus c = new Consensus(argv[0], 0.5);
			c.process();
			c.construct();
			pw.println(c);
			done = true;
		    } catch (DoubleCollision de) {
			System.out.println("Double Collisions found.  Restarting!");
		    }
		}
	    } else {
		System.out.println("Error in opening " + outputfile);
	    }
	} catch (IOException e) {
	    System.out.println("Error in writing " + outputfile);
	    System.exit(-1);
	} finally {
	    pw.close();
	}
    }
}
