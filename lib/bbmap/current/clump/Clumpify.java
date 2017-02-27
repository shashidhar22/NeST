package clump;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import align2.Shared;
import align2.Tools;
import dna.Parser;
import dna.Timer;
import fileIO.ReadWrite;

/**
 * @author Brian Bushnell
 * @date Nov 6, 2015
 *
 */
public class Clumpify {

	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		Clumpify cl=new Clumpify(args);
		cl.process(t);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Clumpify(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		args2=new ArrayList<String>();
		args2.add("in");
		args2.add("out");
		args2.add("groups");
		args2.add("ecco=f");
		args2.add("rename=f");
		args2.add("shortname=f");
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("out") || a.equals("out1")){
				out1=b;
			}else if(a.equals("groups") || a.equals("g") || a.equals("sets")){
				groups=Integer.parseInt(b);
			}else if(a.equals("delete")){
				delete=Tools.parseBoolean(b);
			}else if(a.equals("usetmpdir")){
				useTmpdir=Tools.parseBoolean(b);
			}else if(a.equals("ecco")){
				ecco=Tools.parseBoolean(b);
			}else if(a.equals("compresstemp") || a.equals("ct")){
				if(b!=null & b.equalsIgnoreCase("auto")){forceCompressTemp=forceRawTemp=false;}
				else{
					forceCompressTemp=Tools.parseBoolean(b);
					forceRawTemp=!forceCompressTemp;
				}
			}else if(a.equals("tmpdir")){
				Shared.setTmpdir(b);
			}else if(a.equals("rename") || a.equals("addname")){
				addName=Tools.parseBoolean(b);
			}else if(a.equals("shortname") || a.equals("shortnames")){
				shortName=Tools.parseBoolean(b);
			}else if(a.equals("seed")){
				KmerComparator.setSeed(Long.parseLong(b));
			}else if(a.equals("hashes")){
				KmerComparator.setHashes(Integer.parseInt(b));
			}else{
				args2.add(arg);
			}
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This is called if the program runs with no parameters */
	private void printOptions(){
		throw new RuntimeException("TODO");
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		String[] args=args2.toArray(new String[0]);
		args[2]="groups="+groups;
		if(groups==1){
			args[0]="in="+in1;
			args[1]="out="+out1;
			args[3]="ecco="+ecco;
			args[4]="rename="+addName;
			args[5]="shortname="+shortName;
			KmerSort.main(args);
		}else{
			Random randy=new Random();
			String temp;
			String core="clumpify", path="", extension=".fq";
			if(out1!=null){
				core=ReadWrite.stripToCore(out1);
				path=ReadWrite.getPath(out1);
				extension=ReadWrite.getExtension(out1);
			}
			if(useTmpdir && Shared.tmpdir()!=null){
				temp=Shared.tmpdir()+core+"_temp%_"+Long.toHexString((randy.nextLong()&Long.MAX_VALUE))+extension;
			}else{
				temp=path+core+"_temp%_"+Long.toHexString((randy.nextLong()&Long.MAX_VALUE))+extension;
			}
			
			String comp=ReadWrite.compressionType(temp);
			if(forceCompressTemp && comp==null){temp+=".gz";}
			else if(comp!=null && forceRawTemp){
				temp=temp.substring(0, temp.lastIndexOf('.'));
			}
			args[0]="in="+in1;
			args[1]="out="+temp;
			args[3]="ecco="+ecco;
			args[4]="addname=f";
			args[5]="shortname="+shortName;
			KmerSplit.maxZipLevel=2;
			KmerSplit.main(args);
			
			args[0]="in="+temp;
			args[1]="out="+out1;
			args[3]="ecco=f";
			args[4]="addname="+addName;
			args[5]="shortname=f";
			KmerSort.main(args);
			
			if(delete){
				for(int i=0; i<groups; i++){
					new File(temp.replaceFirst("%", ""+i)).delete();
				}
			}
		}
		t.stop();
		System.err.println("Total time: \t"+t);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int groups=32;
	private boolean ecco=false;
	private boolean addName=false;
	private boolean shortName=false;
	private boolean useTmpdir=true;
	private boolean delete=true;
	private boolean forceCompressTemp=false;
	private boolean forceRawTemp=false;
	
	private String in1=null;
	private String out1=null;
	
	ArrayList<String> args2=new ArrayList<String>();
	private PrintStream outstream=System.err;
	
}
