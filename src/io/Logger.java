package io;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;

public class Logger {

	public enum DebugMode{result, brief, normal, detailed}

	static DebugMode log_level = DebugMode.detailed;
	static OutputStreamWriter logFile;
//	static int counter = 0;
	
	static DateFormat formatter = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss");

	public Logger(){
		try {
			new File("./logs").mkdir();
			logFile = new OutputStreamWriter( new FileOutputStream("logs/log"+System.currentTimeMillis()+".log"),"UTF-16");
		} catch (FileNotFoundException e) {
			System.err.println("No logging available, create a logs folder to record the logs.");
//			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	
	public Logger(String path){
		try {
		//	new File("./logs").mkdir();
			logFile = new OutputStreamWriter( new FileOutputStream(path + "/res"+formatter.format(System.currentTimeMillis())+".txt"),"UTF-16");
		} catch (FileNotFoundException e) {
			System.err.println("No logging available, create a logs folder to record the logs.");
//	ls
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void setLog_level(DebugMode level) {
		log_level = level;
	}

//	public static void log(Object s){
//		try {
//			logFile.write(s.toString().getBytes());
//		} catch (IOException e) {
//			e.printStackTrace();
//		}	
//	}
//	
//	public static void logln(Object s){
//		try {
//			logFile.write((s+"\n").getBytes());
//		} catch (IOException e) {
//			e.printStackTrace();
//		}	
//	}
//	
	public static void log(Object s, DebugMode level){
		if(logFile == null || level.compareTo(log_level) > 0) return;
		try {
			logFile.write(s.toString());
//			logFile.writeUTF(s.toString());
			if(level == DebugMode.result) System.err.print(s);
			logFile.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	public static void log(Object s){
		log(s,DebugMode.normal);	
	}
	public static void logln(Object s){
		logln(s,DebugMode.normal);	
	}
	public static void logln(Object s, DebugMode level){
		if(logFile == null || level.compareTo(log_level) > 0) return;
		try {
			logFile.write((s+"\n"));
//			logFile.writeUTF((s+"\n"));
			if(level == DebugMode.result) System.err.println(s);
			logFile.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	public static void logFunction(Object fnName){
		logFunction(fnName.toString(), DebugMode.normal);
	}
	public static void logFunction(String fnName, DebugMode level){
		if(logFile == null || level.compareTo(log_level) > 0) return;
		try {
			logFile.write(("%----------------------------------------------- "+fnName+" -----------------------------------------------\n"));
			if(level == DebugMode.result) System.err.print("%----------------------------------------------- "+fnName+" -----------------------------------------------\n");
			logFile.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	
}
