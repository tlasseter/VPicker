//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.MVC;

import com.Sts.Actions.Build.StsPassword;
import com.Sts.IO.StsAbstractFile;
import com.Sts.IO.StsFile;
import com.Sts.IO.StsJar;
import com.Sts.UI.Beans.StsJPanel;
import com.Sts.UI.Icons.StsIcon;
import com.Sts.UI.StsMessage;
import com.Sts.UI.StsOkDialog;
import com.Sts.UI.StsPasswordDialog;
import com.Sts.UI.StsYesDialog;
import com.Sts.UI.Toolbars.StsTimeActionToolbar;
import com.Sts.Utilities.StsException;
import com.Sts.Utilities.StsMath;
import com.Sts.Utilities.StsToolkit;
import com.s2s.client.usage.UsageCharge;
import com.s2s.client.usage.UsageException;
import com.s2s.client.usage.UsageManager;

import javax.media.opengl.*;
import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.net.SocketTimeoutException;
import java.util.Calendar;
import java.util.Date;


public class Main
{
    static public String OS;
    static public boolean isDebug = false;
    static boolean debug = false;
    static private boolean suppressMessages = false;
    static public boolean debugPoint = false;
    static public int debugPointRow = 68, debugPointCol = 147;
    static public boolean isWebStart = false;
    static public boolean isJarFile = false;
    static public boolean isDbDebug = false;
    static public boolean isCollabDebug = false;
    static public boolean isDbCmdDebug = false;
    static public boolean isDbIODebug = false;
	static public boolean isDbPosDebug = false;
    static public boolean usageDebug = true;
    static public boolean isGLDebug = false;
    static public boolean isGLTrace = false;
    static public boolean isDrawDebug = false;
	static public boolean isSwingCheck = false;
    static public boolean viewerOnly = false;
    static public double javaVersionNumber = 1.0;
	static public double fovy=50.;
	static public boolean useJPanel = false;
    static public final boolean useNewJavaVersion = false;
    static public final boolean isCollaboration = false;
    static public final boolean isCoreLabs = true;
	public static final String LauchFileFlag = "launchFile";
    static public boolean usageTracking = true;
    static public boolean killSwitch = false;
    static public String usageModule = "com.Sts.MVC.StsMain";
    static public String usageMessage = "None";
    static public String startModule = null;
    static public boolean isModelIntialized = false;
    static private StsModel model = null;
    static public String[] workflowClasses = null;
    static public Boolean[] workflowStatus = null;
    static public String vendorName = "GeoCloudRealTime";
    static public boolean singleWorkflow = false;
    static public StsJar jar = null;

    static public boolean isJarDB = false;
    static public String jarDBFilename = "dbWebStart.jar";
	static public boolean isFileDB = false;
	static public String fileDBPathName;
    static public int timeout = 200;  // Connection timeout to run offline.

    static public boolean usageOK = false;
    static public long usageTime = System.currentTimeMillis(); // milliseconds time since last log

    static public Point screenCenter = new Point(500, 500);

    static public Calendar calendar = Calendar.getInstance();
    static public boolean isJava5 = false;
    static public String[] S2Sargs = null;

	static public GLContext sharedContext;

    //METHOD:  Main method
    public static void main(String[] args)
    {

        // set mainDebug or production flag
        if (args.length > 0) setArgs(args);

        StsToolkit.runWaitOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    createShareableContext();
                    runApplication();
                    System.out.println("launched app");
                }
            }
        );
    }

    private static class Listen implements GLEventListener
    {
        public void init(GLAutoDrawable drawable)
        {

            GL gl = drawable.getGL();
            gl.glClear(GL.GL_COLOR_BUFFER_BIT);
            gl.glFinish();
            Main.sharedContext = drawable.getContext();
            System.err.println("INIT CONTEXT IS: " + drawable.getContext());
        }

        public void reshape(GLAutoDrawable drawable, int x, int y, int width, int height) {}

        public void display(GLAutoDrawable drawable)
        {
            GL gl = drawable.getGL();
            gl.glClear(GL.GL_COLOR_BUFFER_BIT);
        }

        public void displayChanged(GLAutoDrawable drawable, boolean modeChanged, boolean deviceChanged) {}
    }

   private static boolean createShareableContext()
   {
	   /** JBW first try Offscreen rendering buffer */
	   GLPbuffer offScreenBuffer = null;
	   int width = 2;
	   int height = 2;

		   GLCapabilities cap = new GLCapabilities();

		   cap.setRedBits(8);
		   cap.setGreenBits(8);
		   cap.setBlueBits(8);
		   cap.setDoubleBuffered(false);

		   try
		   {
			   offScreenBuffer = GLDrawableFactory.getFactory().createGLPbuffer(cap, null, width, height, null);
			   Listen l = new Listen();
			   offScreenBuffer.addGLEventListener(l);
			   offScreenBuffer.getContext().makeCurrent();
			   offScreenBuffer.getGL().glClear(GL.GL_COLOR_BUFFER_BIT);
			   offScreenBuffer.getGL().glFinish();
			   sharedContext = offScreenBuffer.getContext();
		   }
		   catch(Exception E)
		   {
			   System.err.println("Failed to create offscreen, retry...");
			   E.printStackTrace();
			   //cap.setOffscreenRenderToTextureRectangle(false);
			   try
			   {
				   offScreenBuffer = GLDrawableFactory.getFactory().createGLPbuffer(cap, null, width, height, null);
				   Listen l = new Listen();
				   offScreenBuffer.addGLEventListener(l);
				   offScreenBuffer.getContext().makeCurrent();

			   }
			   catch(Exception E2)
			   {
				   E2.printStackTrace();
				   System.err.println("Failed utterly to create offscreen");
			   }
		   }

		   if(offScreenBuffer == null)
		   {

			   System.err.println("Failed utterly to create offscreen");
			   try
			   {
				   cap.setDoubleBuffered(true);
				   final GLCanvas glc = new GLCanvas(cap, null, null, null);


				   // jbw 3/18 hack to get JOGL to keep the shared context alive;
				   //     the shared context has to be actually realized on the underlying peer
				   glc.setSize(width,height);
				   final JWindow joglFrame = new JWindow();
				   joglFrame.setSize(width,height);
				   joglFrame.add(glc);
				   joglFrame.pack();
				   Listen l = new Listen();
				   glc.addGLEventListener(l);
				   /* blink */
				   joglFrame.setVisible(true);
				   //joglFrame.setVisible(false);

				   return true;
			   }
			   catch(Exception e)
			   {
				   StsException.outputWarningException(null, "createShareableContext", e);
				   e.printStackTrace();
				   return false;
			   }
		   }
		   return true;
   }


    static private void runApplication()
    {
        try
        {
            // get the operating system name
            OS = System.getProperty("os.name");
            System.out.println("OS: " + OS);

            // set the look and feel

            if (OS.startsWith("Windows"))
            {
				try
				{
                    UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
                /*
					for (UIManager.LookAndFeelInfo info : UIManager.getInstalledLookAndFeels())
					{
						if ("Nimbus".equals(info.getName()))
						{
							UIManager.setLookAndFeel(info.getClassName());
							break;
						}
					}
			    */
				}
				catch(Exception e)
				{
					UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.WindowsLookAndFeel());
					System.out.println("Using Windows look-and-feel...");
				}

            }
            else if (OS.equals("Solaris"))
            {
                UIManager.setLookAndFeel(new com.sun.java.swing.plaf.motif.MotifLookAndFeel());
                System.out.println("Using Motif look-and-feel...");
            }
            else if (OS.equals("Irix"))
            {
                UIManager.setLookAndFeel(new com.sun.java.swing.plaf.motif.MotifLookAndFeel());
                System.out.println("Using Motif look-and-feel...");
            }
            else if (OS.equals("Linux"))
            { // Could probably replace all the above with this(?):
                //UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
                //UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.WindowsLookAndFeel());
                UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
                System.out.println("Using Java look-and-feel...");
            }
            else
            {
                // Default: Java look and feel:
                UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
            }
            float test = Float.MAX_VALUE;
            // A hack until we can spend more time on the usage monitoring and operating offline.
            //                    10/29/06               sec/ms  * sec/min * min/hr * hr/day
            long msPerDay = 1000L * 60L * 60L * 24L;
            System.out.println("Current time is: " + System.currentTimeMillis());
            Date killDate = new Date((1304528778913L + (long) (msPerDay * 35L))); //45 days from 05/04/11
            if (killSwitch)
            {
                long killTime = killDate.getTime();
                long currentTime = System.currentTimeMillis();
                if (currentTime > killTime)
                {
                    StsYesDialog.questionValue(new Frame(), "DEMONSTRATION LICENSE HAS EXPIRED\n\n");
                    System.exit(0);
                }
                else
                {
                    if (!StsYesDialog.questionValue(new Frame(), "RUNNING DEMONSTRATION LICENSE\n\n You have " +
                        (float) ((killTime - currentTime) / msPerDay) +
                        " days remaining in demonstration license.\n\n Press OK to continue."))
                        System.exit(0);
                }
            }
            setVersion();

            // Notify the user if the software is configured in anyway other than normal.
            if((viewerOnly) || (isJarDB) || (!vendorName.equalsIgnoreCase("GeoCloudRealTime")))
               modeNotification();

            if (isJarDB)
            {
                /*
                if (isWebStart)
                {
                    jar = StsWebStartJar.constructor(jarDBFilename);
                    if(jar == null)
                    {
                        System.out.println("Failed to find WebStart DB jar file: " + jarDBFilename);
                        model = StsModel.constructSplashMainWindow();
                        return;
                    }
                }

                else
                {
                */
                    if(jarDBFilename == null)
                        jarDBFilename = userSpecifiedJar();
                    else
                    {
                        jar = StsJar.constructor("", jarDBFilename);
                        if(jar == null)
                        {
                            jar = StsJar.getJarInClassPath(jarDBFilename);
                            if (jar == null)
                            {
                                new StsMessage(new Frame(), StsMessage.WARNING, "Unable to find " + jarDBFilename + " in classpath....Please locate and select the file.");
                                jarDBFilename = userSpecifiedJar();
                            }
                            if(jarDBFilename == null)
                                System.exit(0);
                        }
                    }
                //}
                String dbName = Main.jarDBFilename.substring(Main.jarDBFilename.lastIndexOf("\\")+1,Main.jarDBFilename.indexOf(".jar"));
                if(!openJarDatabase("db." + dbName, StsJar.constructor("", jarDBFilename)))
                    System.exit(0);
            }
            else if (isFileDB)
            {
                StsFile file = StsFile.constructor(fileDBPathName);
//               StsAbstractFile abstractFile = jar.getFile(0);
                model = StsModel.constructor(file);
                model.initializeActionStatus();
                model.refreshObjectPanel();
            }
            else
            {
                try
                {
                    model = StsModel.constructSplashMainWindow();
                }
                catch (Exception e)
                {
                    StsException.outputFatalException(Main.class, "StsModel.constructor", e);
                    System.exit(0);
                }
            }


            ToolTipManager toolTipManager = ToolTipManager.sharedInstance();
            toolTipManager.setLightWeightPopupEnabled(false);
            JPopupMenu.setDefaultLightWeightPopupEnabled(false);

            isModelIntialized = true;

        }
        catch (Exception e)
        {
            System.out.println("Exception in Main()\n" + e);
            e.printStackTrace();
        }
        catch (java.lang.OutOfMemoryError e2)
        {
            System.out.println("Out of Memory! in Main()\n" + e2);
            e2.printStackTrace();
        }
    }

    static boolean openJarDatabase(String databaseName, StsJar jarNew)
    {
         jar = jarNew;
         String[] filenames = jar.getFilenames();
         if (filenames == null || filenames.length < 1)
         {
             System.out.println("Failed to find any files in jar file: " + jarDBFilename);
             return false;
         }
         else
         {
             String dbFilename = null;
             int dbFileIndex = -1;
             for(int i=0; i<filenames.length; i++)
             {
                 if(filenames[i].contains(databaseName))
                 {
                     dbFilename = filenames[i];
                     dbFileIndex = i;
                     break;
                  }
             }
             if(dbFileIndex == -1)
             {
                 new StsMessage(new Frame(), StsMessage.WARNING, "Unable to find valid S2S database in selected Java archive (" + jarDBFilename + ").\n\n Exiting.");
                 return false;
             }
             System.out.println("Opening " + dbFilename + " in jar file: " + jarDBFilename);
             StsAbstractFile abstractFile = jar.getFile(dbFileIndex);
             model = StsModel.constructor(abstractFile);
             model.initializeActionStatus();
             model.refreshObjectPanel();
        }
        return true;
    }
    static public String userSpecifiedJar()
    {
         JFileChooser chooser = new JFileChooser(".");
         chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
         chooser.setDialogTitle("Select jar file that contains the database and press the Open button");
         chooser.setApproveButtonText("Open File");
         File selectedFile = null;
         while(selectedFile == null)
         {
            int result = chooser.showOpenDialog(null);
            if(result == JFileChooser.CANCEL_OPTION)
                System.exit(0);
            selectedFile = chooser.getSelectedFile();
         }
         if(selectedFile == null)
         {
            System.out.println("Failed to find jar file in path: " + jarDBFilename);
            return null;
         }
         else
         {
            return selectedFile.getAbsolutePath();
         }
    }
    static public void setVersion()
    {
        String javaVersionString = System.getProperty("java.version");
        System.out.println("Java version: " + javaVersionString);
        String subString = javaVersionString.substring(0, 3);
        double javaVersionNumber = Double.parseDouble(subString);
        isJava5 = javaVersionNumber >= 1.5;
    }

    static public void startUsage()
    {
        boolean isServerOnline = true;
        String reason = null;
        if (usageTracking)
        {
            UsageManager usageManager = null;
            try
            {
                if(usageDebug) System.out.println("startUsage:Get Instance of Usage Manager.");
                usageManager = UsageManager.getInstance(); // (timeout);
                System.out.println("Connection Timeout= " + timeout);
                usageManager.setConnectionTimeout(timeout);
            }
            catch (Exception ex)
            {
                ex.printStackTrace();
//                return;                            // Test this addition......SAJ
            }
            try
            {
                if(usageDebug) System.out.println("startUsage:usageManager.reportToServer");
                usageManager.reportToServer();
            }
            catch (UsageException ue)
            {
                isServerOnline = false;
                reason = "Unable to report to server: Either disconnected from Internet or License problems.\n\n";
            }
            catch(Exception e)
            {
                StsMessageFiles.infoMessage("Connection failed or timed out. Logging locally.");
                e.printStackTrace();
            }
            try
            {
                if(usageDebug) System.out.println("startUsage:usageManager.writetoDisk");
                usageManager.writeToDisk();
            }
            catch (Exception ex)
            {
                reason = "Failed to write license information to local cache.\n\n";
                ex.printStackTrace();
            }

            usageOK = false;
            if (usageManager != null)
            {
                if(usageDebug) System.out.println("startUsage: usageManager.startUsage("+usageModule+","+usageMessage+")");
                usageOK = usageManager.startUsage(usageModule, usageMessage);
            }

//            if((!usageOK) && (isModelIntialized))
            if (!usageOK)
            {
                if (isServerOnline)
                {
                    // License Failed
                    if(reason != null)
                        reason = reason + " " + reasonLicenseFailed(usageManager);
                    else
                        reason = reasonLicenseFailed(usageManager);

                    // Ticket Failed
                    if (reason != null)
                        reason = reason + " " + reasonTicketFailed(usageManager);
                    else
                        reason = reasonTicketFailed(usageManager);

                    if (reason != null)
                        JOptionPane.showMessageDialog(null, "Remote usage monitoring not functioning properly:\n    " + reason);
                }
                else
                    JOptionPane.showMessageDialog(null, "Local usage monitoring not functioning properly:\n    " + reason);
                checkPassword();
            }
        }
    }

    static public boolean checkPassword()
    {
        if(viewerOnly)
        {
            usageTracking = false;
            return true;
        }
        String javaVersionString = System.getProperty("java.version");
        if(javaVersionString.contains("1.5"))
        {
            JOptionPane.showMessageDialog(null, "Password override is not supportted in this version of Java:\n\n   " +
                "Please upgrade to version 1.6+ or Contact S2S (281)684-0576 to\n" +
                "   resolve the license problem or to get a valid password.");
            System.exit(0);
        }
        if(validatePassword(getPasswordFromUser()))
            usageTracking = false;
        else
        {
            JOptionPane.showMessageDialog(null, "Supplied password has either expired or is invalid:\n\n   " +
                "Contact S2S (281)684-0576 to resolve the license problem or\n" +
                "    to get a valid password.");
            System.exit(0);
        }
        return true;
    }

    static public boolean validatePassword(String pword)
    {
        StsPassword passwordTool = null;
        long lastUsed = 0L;
        if ((pword == null) || (pword.length() < 8))
            return false;
        try
        {
            passwordTool = new StsPassword(StsPassword.readSecretKey());
            lastUsed = passwordTool.readLastUsed();
        }
        catch (Exception ex)
        {
            System.out.println("Problem creating password object");
            return false;
        }

        // Check last used
        if (lastUsed > System.currentTimeMillis())
            return false;
        else
            passwordTool.writeLastUsed(StsPassword.readSecretKey());

        // Determine a valid password range
        String dPassword = passwordTool.decrypt(pword);
        if (dPassword == null) return false;

        Date sdate = new Date();
        sdate.setTime(Long.parseLong(dPassword.substring(3)));

        Date eDate = new Date();
        eDate.setTime(sdate.getTime() + 2500000000l); // 1 month

        if ((sdate.before(new Date(System.currentTimeMillis()))) && (eDate.after(new Date(System.currentTimeMillis()))))
            return true;
        else
            return false;
    }

    static public String getPasswordFromUser()
    {
        String password = null;
        try
        {
            StsPasswordDialog dialog = new StsPasswordDialog(null, "Enter Password...", "");
            dialog.setVisible(true);
            password = dialog.getText();
            if (password.length() < 8)
                password = null;
        }
        catch (Exception e)
        {
            StsException.outputException("Main.getPasswordFromUser() failed.", e, StsException.WARNING);
            password = null;
        }
        return password;
    }

    static public String reasonLicenseFailed(UsageManager manager)
    {
        long timeLeft = 0L;
        if(usageDebug)
        {
            System.out.println("Time to license start=" + manager.getTimeToLicenseStart());
            System.out.println("Time since license expired=" + manager.getTimeSinceLicenseExpired());
        }
        if(manager.getTimeToLicenseStart() < 0L) // Giving 10 minute break if clocks are not exact.
        {
            timeLeft = -(manager.getTimeToLicenseStart() / 1000L);
            return new String(timeLeft + " seconds till your license is valid (check that your system clock is accurate), or contact S2S or administrator");
        }
        if (manager.getTimeSinceLicenseExpired() > 0L)
        {
            timeLeft = manager.getTimeSinceLicenseExpired();
            return new String("License expired on " + manager.getLicenseDeathDate() + ", or contact S2S or administrator");
        }
        return null;
    }

    static public String reasonTicketFailed(UsageManager manager)
    {
        long timeLeft = 0L;
        if(usageDebug)
        {
            System.out.println("Time to ticket start=" + manager.getTimeToTicketStart());
            System.out.println("Time since ticket expired=" + manager.getTimeSinceTicketExpired());
            System.out.println("Remaining ticket Usage=" + manager.getRemainingUsage());
        }
        if (manager.getTimeToTicketStart() < 0L)
        {
            timeLeft = -(manager.getTimeToTicketStart() / 1000L);
            return new String(timeLeft + " seconds till your ticket is valid (check that your system clock is accurate), or contact S2S or administrator");
        }

        if (manager.getTimeSinceTicketExpired() > 0L)
            return new String("Ticket expired on " + manager.getTicketDeathDate() + ", re-connect to Internet for new ticket.");

        if (manager.getRemainingUsage() < 0L)
            return new String("Ticket usage limit expired, re-connect to Internet for new ticket.");

        return new String("License or ticket invalid, contact S2S or administrator.");
    }

    static public void logUsageTimer(String module, String message)
    {
//        System.out.println("LogUsageTimer" + module + "----" + message);
        if(model.getProject().isRealtime())
        {
            module = StsTimeActionToolbar.realtimeModule;
            message = "Realtime";
        }

        if (module.equals(usageModule))
            logUsageTimer();
        else
        {
            logUsage();

            usageModule = module;
            usageMessage = message;
            if (usageModule == null)
            {
                usageModule = "com.Sts.WorkflowPlugIn.PlugIns." + model.workflowPlugIn.name;
                usageMessage = "None";
            }
        }
    }

    static public void logUsageTimer()
    {
        // Log Message every five minutes
//        System.out.println("Attempting to log message....");
        if ((System.currentTimeMillis() - usageTime) > 300000)
        {
//            System.out.println("logging message....");
            logUsage();
            usageTime = System.currentTimeMillis();
        }
    }

    static public void logReset()
    {
        if (model == null) return;
        usageModule = "com.Sts.WorkflowPlugIn.PlugIns." + model.workflowPlugIn.name;
        usageMessage = "None";
        usageTime = System.currentTimeMillis();
    }

    static public void setLogModule(String module, String message)
    {
        usageModule = module;
        usageMessage = message;
        if (usageModule == null)
        {
            usageModule = "com.Sts.WorkflowPlugIn.PlugIns." + model.workflowPlugIn.name;
            usageMessage = "None";
        }
    }

    static public void logUsage(String module, String message)
    {

        logUsage();

        usageModule = module;
        usageMessage = message;
        if (usageModule == null)
        {
            usageModule = "com.Sts.WorkflowPlugIn.PlugIns." + model.workflowPlugIn.name;
            usageMessage = "None";
        }
        usageTime = System.currentTimeMillis();
    }

    static public void logUsage(boolean suppress)
    {
        suppressMessages = suppress;
        logUsage();
    }

    static public void logUsage()
    {
        if(model != null)
        {
            if(model.getProject().isRealtime())
            {
                usageModule = StsTimeActionToolbar.realtimeModule;
                usageMessage = "Realtime";
            }
        }
        if (usageModule == null || model == null) return;
//        if(usageModule.equals("com.Sts.WorkflowPlugIn.PlugIns." + model.workflowPlugIn.name)) return;
        if (debug) System.out.println("logUsage: Module: " + Main.usageModule + " Message: " + Main.usageMessage);

        if (!isModelIntialized) return;  // Don't want to log on the initialization of every plugin

        if (usageTracking)
        {
            UsageManager usageManager = null;

            try
            {
                if(usageDebug) System.out.println("logUsage:Get Instance of Usage Manager.");
                usageManager = UsageManager.getInstance(); // (timeout);
                usageManager.setConnectionTimeout(timeout);
            }
            catch (Exception ex)
            {
                ex.printStackTrace();
            }
            try
            {
                if(usageDebug) System.out.println("logUsage:Calling usageManager.enterUsage()");
                UsageCharge charge = usageManager.enterUsage();
            }
            catch (Exception ex)
            {
                ex.printStackTrace();
            }
            try
            {
                if(usageDebug) System.out.println("logUsage:Calling usageManager.reporttoServer()");
                usageManager.reportToServer();
            }
            catch(UsageException ue)
            {
                if (!suppressMessages)
                        StsMessageFiles.logMessage("Unable to report to server: Either disconnected from Internet or License problems.");
            }
            catch(SocketTimeoutException e)
            {
                StsMessageFiles.infoMessage("Connection failed or timed out. Logging locally.");
                e.printStackTrace();
            }
            catch(Exception e)
            {
                StsMessageFiles.infoMessage("Connection exception. Logging locally.");
                e.printStackTrace();
            }
            try
            {
                if(usageDebug) System.out.println("logUsage:Calling usageManager.writeToDisk()");
                usageManager.writeToDisk();
            }
            catch (Exception ex)
            {
                if(!viewerOnly)
                {
                    if (!suppressMessages)
                    {
                        JOptionPane.showMessageDialog(null, "Usage logging to local storage failed...exitting");
                        ex.printStackTrace();
                    }
                    System.exit(0);
                }
            }
            StsMessageFiles.infoMessage("Logging Usage at " + new Date(usageTime) + " for workflow step: " + usageModule);
            // Force usage message to constant until usage monitoring code is altered to not create new plugin when the message changes.
            usageMessage = "Standard Usage";
            boolean usageOK = usageManager.startUsage(usageModule, usageMessage);
            if (!usageOK && !viewerOnly)
            {
                String reason = null;
                reason = reasonLicenseFailed(usageManager);
                if ((reason != null) && (!suppressMessages))
                {
                    JOptionPane.showMessageDialog(null, reason);
                    System.exit(0);
                }
                reason = reasonTicketFailed(usageManager);
                if ((reason != null) && (!suppressMessages))
                {
                    JOptionPane.showMessageDialog(null, reason);
                    System.exit(0);
                }
                System.exit(0);
            }
        }
        else
            StsMessageFiles.infoMessage("Usage Tracking Disabled: Cannot save log at " + new Date(usageTime) + " for workflow step: " + usageModule);
    }

    static public void endUsage()
    {
        if (usageTracking)
        {
            UsageManager usageManager = null;

            if (usageOK)
            {
                try
                {
                    if(usageDebug) System.out.println("endUsage:Get instance of usage manager");
                    usageManager = UsageManager.getInstance(); // (timeout);
                    usageManager.setConnectionTimeout(timeout);
                }
                catch (Exception ex)
                {
                    ex.printStackTrace();
                }
                try
                {
                    if(usageDebug) System.out.println("endUsage:Calling usageManager.enterUsage()");
                    usageManager.enterUsage();
                }
                catch (Exception ex)
                {
                    ex.printStackTrace();
                }
                try
                {
                    if(usageDebug) System.out.println("endUsage:Calling usageManager.reportToServer()");
                    usageManager.reportToServer();
                }
                catch (UsageException ux)
                {
                    StsMessageFiles.logMessage("Unable to report to server: Either disconnected from Internet or License problems.");
                }
                catch(SocketTimeoutException ex)
                {
                    StsMessageFiles.infoMessage("Connection failed or timed out. Logging locally.");
                    ex.printStackTrace();
                }
                catch(Exception ex)
                {
                    StsMessageFiles.logMessage("Unable to report to server: Either disconnected from Internet or License problems.");
                    ex.printStackTrace();
                }
                try
                {
                    if(usageDebug) System.out.println("endUsage:Calling usageManager.writeToDisk()");
                    usageManager.writeToDisk();
                }
                catch (Exception ex)
                {
                    JOptionPane.showMessageDialog(null, "Write to storage failed");
                    ex.printStackTrace();
                }
                usageOK = false;
            }
        }
    }

    static private void setArgs(String[] args)
    {
    	Main.S2Sargs = args;
        int nArgs = args.length;
        int n = 0;
        while (n < nArgs)
        {
            String arg = args[n];
            if (stringContainsString(arg, "workflow"))
            {
                if(argStartsWith(arg, "-"))   // negative sign indicates they can view the contents but not run them.
                {
                    workflowStatus = (Boolean[]) StsMath.arrayAddElement(workflowStatus, new Boolean(false));
                    arg = arg.substring(arg.indexOf("-")+1);
                }
                else
                    workflowStatus = (Boolean[]) StsMath.arrayAddElement(workflowStatus, new Boolean(true));
                workflowClasses = (String[]) StsMath.arrayAddElement(workflowClasses, arg);
            }
            else if (argStartsWith(arg, "mainDebug"))
            {
                isDebug = true;
                System.out.println("Debug is on.");
            }
            else if (argStartsWith(arg, "jardb"))
            {
                isJarDB = true;
                if(args[++n].endsWith(".jar"))
                {
                    jarDBFilename = args[n];
                    System.out.println("Database loading from jar: " + jarDBFilename);
                }
                else
                {
                    jarDBFilename = null;
                    n--;
                    System.out.println("Database loading from user selected jar file.");
                }
            }
			else if (argStartsWith(arg, "filedb"))
			{
				isFileDB = true;
				fileDBPathName = args[++n];
				StsFile file = StsFile.constructor(fileDBPathName);
				if(file.exists() && file.filename.startsWith("db."))
				{
					System.out.println("Database loading from file: " + fileDBPathName);
				}
				else
				{
					System.out.println("Failed to load from file DB: " + fileDBPathName);
					// back up in case user forgot to put argument in and it is part of the next option
					fileDBPathName = null;
					n--;
				}

			}
            else if(stringContainsString(arg, "connecttimeout"))
            {
                timeout = Integer.valueOf(args[++n].trim()).intValue();
                System.out.println("Connection timeout set to " + timeout + " from command line.");
            }
            else if (argStartsWith(arg, "j"))
            {
                isJarFile = true;
                System.out.println("This is a JarFile configuration.");
            }
		/*
            else if (argStartsWith(arg, "proj"))
            {
                isProject = true;
                System.out.println("This is a PROJECT Hack.");
            }
        */
            else if (argStartsWith(arg, "w"))
            {
                isWebStart = true;
                System.out.println("This is a WebStart configuration.");
            }
            else if (argStartsWith(arg, "vendor"))
            {
                vendorName = args[++n];
                System.out.println("Configuring for vendor " + vendorName + ".");
            }
            else if (argStartsWith(arg, "p"))
            {
                if (n + 3 > nArgs)
                {
                    System.out.println("Need row and col values following point mainDebug flag.");
                    continue;
                }
                debugPoint = true;
                debugPointRow = Integer.parseInt(args[++n]);
                debugPointCol = Integer.parseInt(args[++n]);
            }
            else if (argStartsWith(arg, "dbd"))
            {
                isDbDebug = true;
                System.out.println("Database Debug is on.");
            }
			else if (argStartsWith(arg, "fovy"))
			{
				fovy = Double.valueOf(args[++n].trim()).doubleValue();
			}
			else if (argStartsWith(arg, "usejpanel"))
			{
				useJPanel = true;
				System.out.println("Will use lightweight GLJPanel");
			}
            else if (argStartsWith(arg, "dbcmdd"))
            {
                isDbCmdDebug = true;
                System.out.println("DB Cmd Debug is on.");
            }
            else if (argStartsWith(arg, "dbiod"))
            {
                isDbIODebug = true;
                System.out.println("DB IO Debug is on.");
            }
			else if (argStartsWith(arg, "dbposd"))
            {
                isDbPosDebug = true;
                System.out.println("DB Position Debug is on.");
            }
            else if (argStartsWith(arg, "cold"))
            {
                isCollabDebug = true;
                System.out.println("Collaboration Debug is on.");
            }
            else if (argStartsWith(arg, "swingcheck"))
            {
                isSwingCheck = true;
                System.out.println("Swing threadcheck Debug is on.");
            }
            else if (argStartsWith(arg, "gld"))
            {
                isGLDebug = true;
            //    System.setProperty("jogl.debug", "true");
                System.out.println("GL Debug is on.");
            }
            else if (argStartsWith(arg, "glt"))
            {
                isGLTrace = true;
            //    System.setProperty("jogl.debug", "true");
                System.out.println("GL Trace is on.");
            }
            else if (argStartsWith(arg, "dd"))
            {
                isDrawDebug = true;
                System.out.println("Draw Debug is on.");
            }

            else if (argStartsWith(arg, "egasuton"))
            {
                usageTracking = false;
                if(usageDebug)System.out.println("Usage tracking disabled.");
            }
            else if (argStartsWith(arg, "killswitch"))
            {
                usageTracking = false;
                if(usageDebug)System.out.println("Kill switch enabled.....Usage tracking disabled.");
                killSwitch = true;
            }
            else if (argStartsWith(arg, "vieweronly"))
            {
                viewerOnly = true;
                if(usageDebug)System.out.println("Viewer only enabled....No object, workflow, or log panels.");
            }
            else if (argStartsWith(arg, "m"))
            {
                startModule = arg.substring(1, arg.length());
                usageModule = startModule;
            }
            else if (argStartsWith(arg, "sysout"))
            {
                String sysoutFilename = args[++n].replaceFirst("<MS>", new Long(System.currentTimeMillis()).toString());
                PrintStream newOutStream = getPrintStream(sysoutFilename);
                if (newOutStream != null)
                {
                    System.out.println("Redirecting System.out to " + sysoutFilename);
                    System.setOut(newOutStream);
                }
            }
            else if (argStartsWith(arg, "syserrout"))
            {
                String sysoutFilename = args[++n].replaceFirst("<MS>", new Long(System.currentTimeMillis()).toString());
                PrintStream newErrStream = getPrintStream(sysoutFilename);
                if (newErrStream != null)
                {
                    System.out.println("Redirecting System.err to " + sysoutFilename);
                    System.setErr(newErrStream);
                }
            }
            else
                System.out.println("Undefined command line arg: " + args[n]);

            n++;
        }
/*
        if(isWebStart && isJarFile)
        {
            System.out.println("Can't be both WebStart(w) and JarFile(j) in command line arguments: setting to WebStart.");
            isJarFile = false;
        }
*/
    }

    static private PrintStream getPrintStream(String filename)
    {
        boolean autoFlush = true;
        PrintStream aPrintStream = null;
        try
        {
            aPrintStream =
                new PrintStream(
                    new FileOutputStream(filename), autoFlush);
        }
        catch (FileNotFoundException ex)
        {
            return null;
        }
        return aPrintStream;
    }

    static private boolean argStartsWith(String string, String match)
    {
        return string.toUpperCase().toLowerCase().startsWith(match);
    }


    static private boolean stringContainsString(String string, String subString)
    {
        return string.toUpperCase().toLowerCase().indexOf(subString) >= 0;
    }

    static public long getTime(int year, int month, int date, int hour, int minute, int second)
    {
        if (calendar == null) calendar = Calendar.getInstance();
        calendar.set(year, month, date, hour, minute, second);
        return calendar.getTimeInMillis();
    }

    static public void modeNotification()
    {
        StsJPanel panel = new StsJPanel();
        JLabel vendorIcon = new JLabel();
        Icon icon = StsIcon.createIcon(vendorName + "SplashLogo.gif");
		if (icon == null)
		    icon = StsIcon.createIcon("GeoCloudRealTimeSplashLogo.gif");
		panel.setBackground(Color.white);
		panel.setBorder(BorderFactory.createLineBorder(Color.black));
		vendorIcon.setBackground(new Color(252, 250, 252));
		vendorIcon.setBorder(BorderFactory.createEmptyBorder(12, 0, 4, 0));
		vendorIcon.setIcon(icon);
        panel.gbc.anchor = GridBagConstraints.CENTER;
        panel.addEndRow(vendorIcon);

        panel.addEndRow(new JLabel(" "));
        panel.addEndRow(new JLabel("SOFTWARE CONFIGURATION"));

        panel.gbc.anchor = GridBagConstraints.WEST;
        JLabel label1, label2, label3;
        JTextArea area1, area2, area3;
        if(Main.isJarDB)
        {
            label1 = new JLabel("    Archive Loading");
            panel.addEndRow(label1);
            area1 = new JTextArea("      Full versions of the software can create archives of projects. These archives are \n" +
                                  "      delivered to other S2S Software users for review and visual analysis. Archives are \n" +
                                  "      read-only.");
            panel.addEndRow(area1);
            panel.addEndRow(new JLabel(" "));
        }

        if(Main.viewerOnly)
        {
            label2 = new JLabel("    Visual Analysis Only");
            panel.addEndRow(label2);
            area2 = new JTextArea("      Visual analysis limits the application functionality to only 2D, 3D and time-dependent \n" +
                                  "      data analysis of archived and native S2S Projects. If any of the capabilities listed below \n" +
                                  "      are required, contact S2S Systems for an upgrade. Projects are read-only in this mode.");
            panel.addEndRow(area2);
            panel.gbc.anchor = GridBagConstraints.CENTER;
            for(int i=0; i<workflowClasses.length; i++)
                panel.add(new JLabel("      " + workflowClasses[i]));
            panel.addEndRow(new JLabel(" "));
        }

        panel.gbc.anchor = GridBagConstraints.WEST;
        if(!Main.vendorName.equalsIgnoreCase("s2ssystems"))
        {
             label3 = new JLabel("    Partner Release");
             panel.addEndRow(label3);
             area3 = new JTextArea("      Partner releases are configuration specific to S2S Systems Partner companies and may not\n" +
                                   "      represent all available functionality and may include partner proprietary contributions \n" +
                                   "      Contact S2S Systems or the supplier of the software to determine if additional features \n " +
                                   "      of interest are available.");
             panel.addEndRow(area3);
             panel.addEndRow(new JLabel(" "));
        }

         panel.addEndRow(new JLabel(" "));
         StsOkDialog dialog = new StsOkDialog(new Frame(), panel, "S2S Project Archive Load", true);
         return;
    }

}






