import java.util.*;
import java.io.*;

public class readtime {

    private String fn1;
    private String fn2;
    
    public readtime(String s1, String s2) {

	fn1 = s1;
	fn2 = s2;

    }

    private void process() throws Exception {

	BufferedReader s1 = new BufferedReader(new FileReader(fn1));
	PrintWriter s2 = new PrintWriter(new FileOutputStream(fn2));
	
	String line1 = null;
	
	while((line1=s1.readLine())!=null)
	    {
		StringTokenizer st1 = new StringTokenizer(line1);
		
		String docid1 = "";
		if (st1.hasMoreElements())
		    docid1 = st1.nextToken();
		
		if (docid1.equals("elapsed_time"))
		    {
			line1=s1.readLine();
			line1=s1.readLine();
			st1 = new StringTokenizer(line1);
			float f1 = Float.parseFloat(st1.nextToken());
			line1=s1.readLine();
			line1=s1.readLine();
			line1=s1.readLine();
			line1=s1.readLine();
			line1=s1.readLine();
			st1 = new StringTokenizer(line1);
			float f2 = Float.parseFloat(st1.nextToken());
			System.out.println("train: "+f1+" test: "+f2);
			s2.println(f1+" "+f2);
		    }
		
	   }
	s1.close();
	s2.close();

    }

    public static void main(String[] args) throws Exception { 

	if (args.length != 2)
	{
	    System.err.println("Usage: java readtime [time_file] [outfile]");
	    System.exit(1);
	}

	readtime p = new readtime(args[0], args[1]);
	p.process();

    }

}
