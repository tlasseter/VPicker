//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Help;

/**
 This class is used to manage help systems using an HTML
 browser and/or JavaHelp
 */

import com.Sts.MVC.*;
import com.Sts.Utilities.*;
import com.Sts.UI.StsButton;

import javax.help.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionListener;
import java.net.*;

public class HelpManager
{
    static private HelpManager helpManager = null;   // singleton pattern
    static private final String GENERAL_HELPSET_FILE = "General/S2SGeneral.hs";
    static public final byte GENERAL = 0;
    static private final String STEPS_HELPSET_FILE = "Workflows/S2SWorkflows.hs";
    static public final byte STEPS = 1;
    static public final byte WIZARD = 1;
    static private final String subdirectory = "";
//    static private final String helpJarname = "S2SHelp.jar";

    private HelpSet helpSet;
    private HelpBroker helpBroker;
    private HelpSet generalHelpSet;
    private HelpBroker generalHelpBroker;
    private HelpSet stepsHelpSet;
    private HelpBroker stepsHelpBroker;
    private ActionListener stepsActionListener;
    private HelpSet wizardHelpSet;
    private HelpBroker wizardHelpBroker;

    static private final boolean debug = false;

    /**
     * get the help manager
     */
    static public HelpManager getHelpManager(Component parentComponent) throws HelpSetException, MalformedURLException
    {
        if (helpManager == null) helpManager = new HelpManager(parentComponent);
        return helpManager;
    }

    static public void deleteManager() { helpManager = null; }

    /** constructor can't be called directly.  It is executed by getHelpManager() */
    protected HelpManager(Component parentComponent) throws HelpSetException, MalformedURLException
    {
        if(!SwingUtilities.isEventDispatchThread())
            System.out.println("help not on swing thread");
        initJavaHelp(parentComponent);
        helpManager = this;
    }

    public boolean addHelpSet(String helpSetName)
    {
        try
        {
            ClassLoader classLoader = getClass().getClassLoader();
            URL hsUrl = HelpSet.findHelpSet(classLoader, subdirectory + helpSetName);
            if (hsUrl == null)
                return false;
            //stepsHelpSet.add(hs);
            wizardHelpSet = new HelpSet(classLoader, hsUrl);
            wizardHelpBroker = wizardHelpSet.createHelpBroker();
        }
        catch (Exception ex)
        {
            System.err.println("Error loading helpset at URL:" + helpSetName + " ex=" + ex);
            ex.printStackTrace();
        }
        return true;
    }

    public boolean initWorkflowHelpSets(JMenuItem menuItem)
    {
        URL stepsHelpURL;
        try
        {
            menuItem.removeActionListener(stepsActionListener);
            stepsHelpBroker = null;
            stepsHelpSet = null;
            stepsActionListener = null;

            ClassLoader classLoader = getClass().getClassLoader();
            stepsHelpURL = HelpSet.findHelpSet(classLoader, subdirectory + STEPS_HELPSET_FILE);
            if (stepsHelpURL == null)
            {
                StsException.systemError("Unable to load workflow steps help files: S2SNewHelp.jar unavailable.");
                return false;
            }
            if (debug) System.out.println("Help url: " + stepsHelpURL.toString());
            stepsHelpSet = new HelpSet(classLoader, stepsHelpURL);
            stepsHelpBroker = stepsHelpSet.createHelpBroker();
            if (helpManager.stepsHelpBroker == null)
                return false;
            stepsActionListener = new CSH.DisplayHelpFromSource(stepsHelpBroker);
            menuItem.addActionListener(stepsActionListener);
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
        return true;
    }

    public boolean addWorkflowHelpSet(String helpSetName)
    {
        try
        {
            if (stepsHelpSet == null)
                return false;
            ClassLoader classLoader = getClass().getClassLoader();
            URL hsUrl = HelpSet.findHelpSet(classLoader, subdirectory + helpSetName);
            if (hsUrl == null)
                return false;
            HelpSet hs = new HelpSet(classLoader, hsUrl);
            if (hs != null)
            {
                stepsHelpSet.remove(hs);
                stepsHelpSet.add(hs);
            }
        }
        catch (Exception ex)
        {
            System.err.println("Error loading helpset at URL:" + helpSetName + " ex=" + ex);
            ex.printStackTrace();
        }
        return true;
    }

    // classInitialize the help system

    private void initJavaHelp(Component parentComponent) throws HelpSetException, MalformedURLException
    {
        URL jarURL, generalHelpURL, stepsHelpURL;
        ClassLoader classLoader = null;

        try
        {
            if (Main.isWebStart)
            {
                URL codebase = JNLPUtilities.getCodeBase();
                String codebaseName = codebase.toString();

                // If help jar is already cached, then use it.
                // if not, ask the user whether s/he wishes to download jar.
                // If not, access directly over the web.
                // If working from a jar, use classloader to access jar entries;
                // otherwise classloader is null (accessing over the web).

                classLoader = getClass().getClassLoader();
                // This section does not work because of the encrypted paths in name
/*
                jarURL = new URL(codebase, helpJarname);
                if(JNLPUtilities.isResourceCached(jarURL))
                {
                    classLoader = getClass().getClassLoader();
                    if(mainDebug) System.out.println("Help jar cached. url: " + jarURL.toString());
                }
                else
                {
                    boolean loadHelpJars = StsMessage.questionValue(parentComponent,
                            "Do you wish to download S2SHelp.jar file?\n" +
                            "If no, you will access them over the web.");
                    if(loadHelpJars)
                    {
                        jarURL = new URL("jar:" + codebaseName + helpJarname + "!/");
                        //if(mainDebug)
                        System.out.println("Help jar not cached.  Downloading url: " + jarURL.toString());
                        JNLPUtilities.downloadJarURL(helpJarname);
                        classLoader = getClass().getClassLoader();
                    }
                    else
                    {
                        classLoader = null;
                        if(mainDebug) System.out.println("Help jar not cached. Not downloading. Will download on demand from url: " + jarURL.toString());
                    }
                }
*/
//                 URLClassLoader urlClassLoader = new URLClassLoader(new URL[] { jarURL } );
//                 helpURL = HelpSet.findHelpSet(urlClassLoader, subdirectory + HELPSET_FILE);
//                 helpSet = new HelpSet(urlClassLoader, helpURL);
            }
            if (debug) System.out.println("General Help file: " + subdirectory + GENERAL_HELPSET_FILE);
            generalHelpURL = HelpSet.findHelpSet(classLoader, subdirectory + GENERAL_HELPSET_FILE);
            if (generalHelpURL == null)
            {
                StsException.systemError("Unable to load general help files: S2SHelp.jar unavailable.");
                return;
            }
            if (debug) System.out.println("Help url: " + generalHelpURL.toString());
            generalHelpSet = new HelpSet(classLoader, generalHelpURL);
            generalHelpBroker = generalHelpSet.createHelpBroker();
            /*
            if(debug) System.out.println("Steps Help file: " + subdirectory + GENERAL_HELPSET_FILE);
            stepsHelpURL = HelpSet.findHelpSet(classLoader, subdirectory + STEPS_HELPSET_FILE);
            if(stepsHelpURL == null)
            {
                StsException.systemError("Unable to load workflow steps help files: S2SHelp.jar unavailable.");
                return;
            }
            if(debug) System.out.println("Help url: " + stepsHelpURL.toString());
            stepsHelpSet = new HelpSet(classLoader, stepsHelpURL);
            stepsHelpBroker = stepsHelpSet.createHelpBroker();
            */
        }
        catch (Exception e)
        {
            StsException.outputException("HelpManager.initJavaHelp() failed.",
                    e, StsException.WARNING);
        }
    }

    // classInitialize the help system

    private void initHelp(Component parentComponent, String helpFile, String helpJar) throws HelpSetException, MalformedURLException
    {
        URL helpURL, jarURL;
        ClassLoader classLoader = null;

        try
        {
            if (Main.isWebStart)
            {
                URL codebase = JNLPUtilities.getCodeBase();
                String codebaseName = codebase.toString();
                if (debug) System.out.println("codebaseName=" + codebaseName);

                // If help jar is already cached, then use it.
                // if not, ask the user whether s/he wishes to download jar.
                // If not, access directly over the web.
                // If working from a jar, use classloader to access jar entries;
                // otherwise classloader is null (accessing over the web).

                jarURL = new URL(codebase, helpJar);
                if (JNLPUtilities.isResourceCached(jarURL))
                {
                    classLoader = getClass().getClassLoader();
                    if (debug) System.out.println("Help jar cached. url: " + jarURL.toString());
                }
/*                else
                {
                    boolean loadHelpJars = StsMessage.questionValue(parentComponent,
                            "Local manual not available, Do you wish to download " + helpJar + " file?\n");
                    if(loadHelpJars)
                    {
                        jarURL = new URL("jar:" + codebaseName + helpJarname + "!/");
                        if(mainDebug) System.out.println("Help jar not cached.  Downloading url: " + jarURL.toString());
                        JNLPUtilities.downloadJarURL(helpJar);
                        classLoader = getClass().getClassLoader();
                    }
                    else
                    {
                        classLoader = null;
                        if(mainDebug) System.out.println("Help jar not cached. Not downloading. Will download on demand from url: " + jarURL.toString());
                    }
                }
 */
            }
            if (debug) System.out.println("Help file: " + subdirectory + helpFile);

            helpURL = HelpSet.findHelpSet(classLoader, subdirectory + helpFile);
            if (helpURL == null)
            {
                StsException.systemError("Unable to load solutions help files: " + helpFile + " unavailable.");
                return;
            }
            if (debug) System.out.println("Help url: " + helpURL.toString());
            helpSet = new HelpSet(classLoader, helpURL);
            if (helpSet == null)
                if (debug) System.out.println("helpSet is null.");

            helpBroker = helpSet.createHelpBroker();
        }
        catch (Exception e)
        {
            StsException.outputException("HelpManager.initJavaHelp() failed.",
                    e, StsException.WARNING);
        }
    }

/*
    public class HelpSetResourceAnchor
    {
       // do nothing; add it to the jar with all your helpsets
    }

    public class HelloHelp extends JFrame
    {
      JHelp   _browser;
      HelpSet _helpset;

      public HelloHelp( String topic_id )
      {
       try
       {
         // use resource anchor trick to reference helpset
         ClassLoader cl = HelpSetResourceAnchor.class.getClassLoader();
         URL url = cl.getResource( "hello.hs" );
         _helpset = new HelpSet( cl, url );

         _browser = new JHelp( _helpset );
         _browser.setNavigatorDisplayed(true);

         getContentPane().setLayout(new GridLayout());
         getContentPane().add( _browser);

         setSize(500, 600);
         centerComponentOnScreen(0, 0);
         setTitle( "Hello Help");

         showHelp( topic_id );

         setVisible(true);
      }
      catch (Exception ex)
      {
         System.out.println("*** error: " + ex.toString() );
      }
    }

    public void showHelp( String topic_id )
    {
      try
      {

       if( topic_id.equals( "" ) )
          topic_id = _helpset.getHomeID();

       _browser.setCurrentID( topic_id );
      }
      catch (Exception ex)
      {
         System.out.println("*** error: " + ex.toString() );
      }
    }
*/

    /**
     * set action listener for source help
     */
    static public boolean setSourceHelpActionListener(StsButton button, byte type, Component component)
    {
        if (button == null) return false;
        try
        {
            getHelpManager(component);
        }
        catch (Exception e)
        {
            return false;
        }

        if (helpManager.wizardHelpBroker == null)
            return false;
        button.addActionListener(new CSH.DisplayHelpFromSource(helpManager.wizardHelpBroker));
        return true;
    }

    /**
     * set action listener for source help
     */
    static public boolean setSourceHelpActionListener(JMenuItem menuItem, byte type, Component component)
    {
        if (menuItem == null) return false;
        try
        {
            getHelpManager(component);
        }
        catch (Exception e)
        {
            return false;
        }
        if (type == GENERAL)
        {
            if (helpManager.generalHelpBroker == null)
                return false;
            menuItem.addActionListener(new CSH.DisplayHelpFromSource(helpManager.generalHelpBroker));
        }
        /*
        else if(type == STEPS)
        {
            if (helpManager.stepsHelpBroker == null)
                return false;
            menuItem.addActionListener(new CSH.DisplayHelpFromSource(helpManager.stepsHelpBroker));
        }
        */
        return true;
    }

    /**
     * set action listener for source help
     */
    static public boolean setSourceHelpActionListener(JMenuItem menuItem, String helpFile, String jarFile, Component component)
    {
        if (menuItem == null)
            return false;
        try
        {
            if (helpManager == null)
                getHelpManager(component);
            helpManager.initHelp(component, helpFile, jarFile);
        }
        catch (Exception e)
        {
            return false;
        }

        if (helpManager.helpBroker == null)
        {
            if (debug) System.out.println("helpManager.helpBroker is null");
            return false;
        }
        menuItem.addActionListener(new CSH.DisplayHelpFromSource(helpManager.helpBroker));
        return true;
    }

    /**
     * set action listener for source help
     */
    static public boolean getRemoteHelp(JMenuItem menuItem, String helpFile, String jarFile, Component component)
    {
        if (menuItem == null)
            return false;
        try
        {
            getHelpManager(component);
            helpManager.initHelp(component, helpFile, jarFile);
        }
        catch (Exception e)
        {
            return false;
        }

        if (helpManager.helpBroker == null)
        {
            if (debug) System.out.println("helpManager.helpBroker is null");
            return false;
        }
        new CSH.DisplayHelpFromSource(helpManager.helpBroker);
        return true;
    }

    /**
     * set action listener for tracking help
     */
    static public boolean setTrackingHelpActionListener(JButton b, byte type, Component component)
    {
        if (b == null) return false;
        try
        {
            getHelpManager(component);
        }
        catch (Exception e)
        {
            return false;
        }
        if (type == GENERAL)
        {
            if (helpManager.generalHelpBroker == null) return false;
            b.addActionListener(new CSH.DisplayHelpAfterTracking(helpManager.generalHelpBroker));
        }
        /*
        else if(type == STEPS)
        {
            if(helpManager.stepsHelpBroker == null) return false;
            b.addActionListener(new CSH.DisplayHelpAfterTracking(helpManager.stepsHelpBroker));
        }
        */
        return true;
    }

    /**
     * set help target button.  Used only in JavaHelp.
     */
    static public boolean setButtonHelp(JButton button, byte type, String target, Component component)
    {
        if (target == null) return true;
        try
        {
            getHelpManager(component);
        }
        catch (Exception e)
        {
            return false;
        }
        if (type == GENERAL)
        {
            if (helpManager.generalHelpBroker == null) return false;
            try
            {
                helpManager.generalHelpBroker.enableHelpOnButton(button, target, null);
            }
            catch (Exception e)
            {
                return false;
            }
        }
        /*
        else if(type == STEPS)
        {
            if(helpManager.stepsHelpBroker == null) return false;
            try { helpManager.stepsHelpBroker.enableHelpOnButton(button, target, null); }
            catch (Exception e) { return false; }
        }
        */
        return true;
    }

    /**
     * set context-sensitive help for UI component
     */
    static public boolean setContextSensitiveHelp(Component comp, byte type, String target, Component component)
    {
        if (comp == null) return false;
        if (target == null) return true;
        try
        {
            getHelpManager(component);
        }
        catch (Exception e)
        {
            return false;
        }
        if (type == GENERAL)
        {
            if (helpManager.generalHelpBroker == null) return false;
            helpManager.generalHelpBroker.enableHelp(comp, target, null);
        }
        /*
        else if(type == STEPS)
        {
            if(helpManager.stepsHelpBroker == null) return false;
            helpManager.stepsHelpBroker.enableHelp(comp, target, null);
        }
        */
        //CSH.setHelpIDString(comp, target);
        return true;
    }

    public static void main(String[] args)
    {
        try
        {
            HelpManager helpManager = getHelpManager(null);
            JHelp browser = new JHelp(helpManager.generalHelpSet);
            browser.setNavigatorDisplayed(true);
//            browser.setCurrentID("");

            JFrame frame = new JFrame("Help Test");
            frame.getContentPane().setLayout(new GridLayout());
            frame.getContentPane().add(browser);
            frame.setSize(500, 600);
            frame.setLocation(0, 0);
            frame.setVisible(true);
        }

        catch (Exception e)
        {
            System.out.println("Unable to classInitialize JavaHelp");
        }
    }
}
