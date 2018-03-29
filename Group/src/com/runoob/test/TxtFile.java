package com.runoob.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class TxtFile {
	private static   String savefile = "D:\\Programme\\Python\\first\\test.txt";  

	public static void saveAsFileWriter(String content) {  
	  
	 FileWriter fw = null;  
	 try {  
	 File f = new File(savefile);
	  fw = new FileWriter(f,true);  
	 } catch (IOException ex) {  
	  ex.printStackTrace();  
	 } 
	 PrintWriter pw = new PrintWriter(fw);
	 pw.print(content);
	 pw.flush();
	 try {
		 fw.flush();
		 pw.close();
		 fw.close();
	 }catch(IOException e) {
		 e.printStackTrace();
	 }
	 }
}
