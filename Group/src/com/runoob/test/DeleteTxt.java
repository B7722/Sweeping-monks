package com.runoob.test;

import java.io.File;

public class DeleteTxt {

	private static   String pfile = "D:\\Programme\\Python\\first\\test.txt";  
	public static void delete() {  
		File f = new File(pfile);
		if(f.exists()) {
			f.delete();
		}
	}
}
