
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;

import com.Sts.MVC.*;
import com.Sts.UI.Sounds.*;

import javax.swing.*;

public class StsException extends Exception
{
    int level;  // exception level: WARNING or FATAL
	public static final int WARNING = 1;
    public static final int FATAL = 2;
    public static final int DEBUG = 3;

    static StsStatusArea statusArea = null;

    static int maxTraceDumps = 10; /** @param maxTraceDumps max number of full traces written */
    static int traceDumps = 0;     /** @param traceDumps current number of traces written */
    static int maxNErrorMessages = 1000; /** @param maxNErrorMessages max number of errors allowed */
    static int nErrorMessages = 0; /** @param nErrorMessages current number of errors   */
    static boolean debug = false;

	static boolean useSounds = false;

    static public void setDebug(boolean staticDebug) { debug = staticDebug; }

	public StsException(int level)
	{
    	super();
    	this.level = level;
	}

    public StsException(int level, String s)
    {
    	super(s);
    	this.level = level;
    }

    public StsException(int level, String id, String s)
    {
    	super( new String(id + ". " + s) );
    	this.level = level;
    }

    public StsException(int level, String id, String s, int value)
    {
    	super( new String(id + ". " + s + " " + String.valueOf(value)) );
    	this.level = level;
    }

    public StsException(int level, String id, String s1, int i1, String s2, int i2)
    {
    	super(new String(id + ". " + s1 + " " + String.valueOf(i1) + " " + s2 + " " + String.valueOf(i2)) );
    	this.level = level;
    }

    public StsException(int level, String id, String s1, int i1, String s2, int i2, String s3, int i3)
    {
    	super(new String(id + ". " + s1 + " " + String.valueOf(i1) + " " +
        	  s2 + " " + String.valueOf(i2) + " " + s3 + " " + String.valueOf(i3)) );
    	this.level = level;
    }

// Methods from here down to mark should be phased out for ones above
// in which we have added the exception level.

	public StsException()
	{
    	super();
	}

    public StsException(String s)
    {
    	super();
    }

    public StsException(String id, String s)
    {
    	super( new String(id + ". " + s) );
    }

    public StsException(String id, String s, int value)
    {
    	super( new String(id + ". " + s + " " + String.valueOf(value)) );
    }

    public StsException(String id, String s1, int i1, String s2, int i2)
    {
    	super(new String(id + ". " + s1 + " " + String.valueOf(i1) + " " + s2 + " " + String.valueOf(i2)) );
    }

    public StsException(String id, String s1, int i1, String s2, int i2, String s3, int i3)
    {
    	super(new String(id + ". " + s1 + " " + String.valueOf(i1) + " " +
        	  s2 + " " + String.valueOf(i2) + " " + s3 + " " + String.valueOf(i3)) );
    }

// Accessors

    public int getLevel() { return level; }

// Mark

    public static void outputException(Exception e, int level)
    {
        outputException(new String(""), e, level);
    }

	static public void outputWarningException(Object object, String methodName, String message, Exception e)
	{
		String error = "Exception: " + StsToolkit.getSimpleClassname(object) + "." + methodName + " failed. " + message;
		outputException(error, e, WARNING);
	}

	static public void outputWarningException(Object object, String methodName, Exception e)
	{
		String error = "Exception: " + StsToolkit.getSimpleClassname(object) + "." + methodName + " failed. ";
		outputException(error, e, WARNING);
	}

	static public void outputWarningException(Class c, String methodName, Exception e)
	{
		String error = "Exception: " + StsToolkit.getSimpleClassname(c) + "." + methodName + " failed. ";
		outputException(error, e, WARNING);
	}

	static public void outputWarningException(Class c, String methodName, String message, Exception e)
	{
		String error = "Exception: " + StsToolkit.getSimpleClassname(c) + "." + methodName + " failed. " + message;
		outputException(error, e, WARNING);
	}

    static public void outputFatalException(Object object, String methodName, String message, Exception e)
	{
		String error = "Exception: " + StsToolkit.getSimpleClassname(object) + "." + methodName + " failed. " + message;
		outputException(error, e, FATAL);
	}

    static public void outputFatalException(Object object, String methodName, Exception e)
	{
		String error = "Exception: " + StsToolkit.getSimpleClassname(object) + "." + methodName + " failed. ";
		outputException(error, e, FATAL);
	}

    static public void outputFatalException(Class c, String methodName, Exception e)
	{
		String error = "Exception: " + StsToolkit.getSimpleClassname(c) + "." + methodName + " failed. ";
		outputException(error, e, FATAL);
	}

	static public void outputWarningException(Object object, String methodName, String message, Error e)
	{
		String errorMessage = "Exception: " + StsToolkit.getSimpleClassname(object) + "." + methodName + " failed. " + message;
		outputException(errorMessage, e, WARNING);
	}

    public static void outputException(final String error, final Exception e, int level)
    {
        // Not doing anything with error string currently

        if(nErrorMessages >= maxNErrorMessages) return;
        nErrorMessages++;

        System.err.println(error);
        StsMessageFiles.errorMessage(error + " " + e.getMessage());
        printStackTrace(e);
        if( level == WARNING)
        {
            StsMessageFiles.logMessage("WARNINGS generated.  See Error dialog.");
            if(useSounds) StsSound.play(StsSound.CAMERA_CLICK);
        }
     	else if(level == FATAL)
        {
 //           sound.setType(StsSound.TOILET);
			StsToolkit.runLaterOnEventThread(new Runnable()
			{
				public void run()
				{

					int answer = JOptionPane.showConfirmDialog(
						null, "Fatal errors!\n"+error+"\n"+e+"\n\nOK to end or Cancel to continue.",
						"System Error", JOptionPane.OK_CANCEL_OPTION,
						JOptionPane.ERROR_MESSAGE);
					if(answer == JOptionPane.OK_OPTION)
					{
						StsMessageFiles.logMessage("Fatal Error.  Exiting...");
						StsSound.play(StsSound.IMPLOSION);
						System.exit(0);
					}
				}
			});
        }
        else if(level == DEBUG && debug)
            StsMessageFiles.logMessage("DEBUG WARNINGS generated.  See Error dialog.");
    }

    public static void outputException(final String errorMessage, final Error e, int level)
    {
        // Not doing anything with errorMessage string currently

        if(nErrorMessages >= maxNErrorMessages) return;
        nErrorMessages++;

        System.err.println(errorMessage);
        StsMessageFiles.errorMessage(errorMessage + " " + e.getMessage());
        e.printStackTrace();
        if(level == WARNING)
        {
            StsMessageFiles.logMessage("WARNINGS generated.  See Error dialog.");
            if(useSounds) StsSound.play(StsSound.CAMERA_CLICK);
        }
     	else if(level == FATAL)
        {
 //           sound.setType(StsSound.TOILET);
			StsToolkit.runLaterOnEventThread(new Runnable()
			{
				public void run()
				{

					int answer = JOptionPane.showConfirmDialog(
						null, "Fatal errors!\n"+ errorMessage +"\n"+e+"\n\nOK to end or Cancel to continue.",
						"System Error", JOptionPane.OK_CANCEL_OPTION,
						JOptionPane.ERROR_MESSAGE);
					if(answer == JOptionPane.OK_OPTION)
					{
						StsMessageFiles.logMessage("Fatal Error.  Exiting...");
						StsSound.play(StsSound.IMPLOSION);
						System.exit(0);
					}
				}
			});
        }
        else if(level == DEBUG && debug)
            StsMessageFiles.logMessage("DEBUG WARNINGS generated.  See Error dialog.");
    }


    static public void printStackTrace(Exception e)
    {
        if(traceDumps < maxTraceDumps)
        {
            e.printStackTrace();
            traceDumps++;
        }
    }
/*
    static private void getStatusArea()
    {
        StsModel currentModel = StsSerializable.getCurrentModel();
        if(currentModel == null) return;
        statusArea = currentModel.getStatusArea();
    }
*/
    static public void systemError(String string)
    {
        if(nErrorMessages >= maxNErrorMessages) return;
        nErrorMessages++;

        String msg = "system error: " + string;
		System.err.println(msg);
        StsMessageFiles.errorMessage(msg);
    }

	static public void systemError(Object object, String methodName)
	{
		if(nErrorMessages >= maxNErrorMessages) return;
		nErrorMessages++;

		String msg = "system error: " + StsToolkit.getSimpleClassname(object) + "." + methodName + " failed.";
        System.err.println(msg);
        StsMessageFiles.errorMessage(msg);
	}

	static public void systemError(Object object, String methodName, String message)
	{
		if(nErrorMessages >= maxNErrorMessages) return;
		nErrorMessages++;

		String msg = "system error: " + StsToolkit.getSimpleClassname(object) + "." + methodName + " error. " + message;
        System.err.println(msg);
        StsMessageFiles.errorMessage(msg);
	}

    static public void systemError(Class c, String methodName)
	{
        if(nErrorMessages >= maxNErrorMessages) return;
		nErrorMessages++;

		String msg = "system error: " + StsToolkit.getSimpleClassname(c) + "." + methodName + " failed.";
        System.err.println(msg);
        StsMessageFiles.errorMessage(msg);
	}

    static public void systemError(Class c, String methodName, String message)
	{
        if(nErrorMessages >= maxNErrorMessages) return;
		nErrorMessages++;

		String msg = "system error: " + StsToolkit.getSimpleClassname(c) + "." + methodName + " failed. " + message;
        System.err.println(msg);
        StsMessageFiles.errorMessage(msg);
	}

	static public void notImplemented(Object object, String methodName)
	{
		if(nErrorMessages >= maxNErrorMessages) return;
		nErrorMessages++;

		String msg = "system error: " + methodName + " not implemented in class " + StsToolkit.getSimpleClassname(object);
		if(!StsMessageFiles.errorMessage(msg))
			System.err.println(msg);
	}

	static public void notImplemented(Class c, String methodName)
	{
		if(nErrorMessages >= maxNErrorMessages) return;
		nErrorMessages++;

		String msg = "system error: " + methodName + " not implemented in class " + StsToolkit.getSimpleClassname(c);
		if(!StsMessageFiles.errorMessage(msg))
			System.err.println(msg);
	}

    static public void systemDebug(String string)
    {
//        if(nErrorMessages >= maxNErrorMessages) return;
//        nErrorMessages++;

        String msg = "system mainDebug: " + string;
        System.out.println(msg);
        StsMessageFiles.errorMessage(msg);
    }

	static public void systemDebug(Object object, String methodName, String message)
	{
//		if(nErrorMessages >= maxNErrorMessages) return;
//		nErrorMessages++;

		String msg = "system debug: " + StsToolkit.getSimpleClassname(object) + "." + methodName + "() " + message;
	    StsMessageFiles.errorMessage(msg);
        System.out.println(msg);
    }

    static public void systemDebug(Class c, String methodName, String message)
     {
         String msg = "system debug: " + StsToolkit.getSimpleClassname(c) + "." + methodName + " () " + message;
         System.out.println(msg);
         StsMessageFiles.errorMessage(msg);
     }

	static public void systemDebug(Object object, String methodName)
	{
//		if(nErrorMessages >= maxNErrorMessages) return;
//		nErrorMessages++;

		String msg = "system debug: " + StsToolkit.getSimpleClassname(object) + "." + methodName + "() ";
        System.out.println(msg);
    }

    public static void outputException(StsException e)
    {
        outputException(e, e.getLevel());
    }
/*
    public void output()
    {
        if(nErrorMessages >= maxNErrorMessages) return;
        nErrorMessages++;

        String msg = "StsException: " + getMessage();
   		System.err.println(msg);
        errorMessage(msg);

     	if(level == StsException.FATAL)
        {
       		StsException.printStackTrace(this);
            logMessage("Fatal Error.  Exiting...");
			errorMessage("Fatal Error.  Exiting...");
        	System.exit(0);
		}
    }
*/
// These methods should be phased out (down to mark) for outputException above

    public void print()
    {
        String msg = "StsException: " + getMessage();
   		System.err.println(msg);
        errorMessage(msg);
    }

	public void errorMessage(String msg)
	{
		StsMessageFiles.errorMessage(msg);
	}

	public void logMessage(String msg)
	{
		StsMessageFiles.logMessage(msg);
	}

    public void warning()
    {
    	print();
	}

    public void fatalError()
    {
    	print();
        printStackTrace(this);
        errorMessage("Fatal Error.  Exiting...");
        logMessage("Fatal Error.  Exiting...");
        System.exit(0);
	}

// Mark

    public static void main(String argv[])
    {
    	int i;

        try
        {
        	i = Integer.parseInt(argv[0]);
        }
        catch (ArrayIndexOutOfBoundsException e)
        {
        	System.err.println("Must specify an argument");
            return;
        }
        catch (NumberFormatException e)
        {
        	System.err.println("Must specify an integer argument.");
            return;
        }

        try
        {
        	methodA(i);
        }
        catch (StsException e)
        {
        	e.print();
        }
        finally
        {
        	System.err.println("Terminating test!");
        }
    }

    public static void methodA(int i) throws StsException
    {
    	switch(i)
        {
        	case 0:
            	throw new StsException("methodA");
            case 1:
            	throw new StsException("methodA", "String and value: ", 1);
            case 2:
            	throw new StsException("methodA", "String one and value:", 1, "string two and value:", 2);
           	case 3:
            {
            	try
                {
                	float[] f = new float[-i];
                }
                catch (Exception e)
                {
            		throw new StsException("methodA", "Float alloc error for size :", i);
                }
            }
            default:
            	throw new StsException("methodA", "Error test not specified for this integer:", i);
        }
    }
}

