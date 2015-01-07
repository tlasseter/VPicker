package com.Sts.DBTypes;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */

import com.Sts.DB.DBCommand.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.magician.fonts.*;

import javax.media.opengl.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;

public class StsAncillaryData extends StsMainTimeObject implements StsTreeObjectI, StsViewSelectable
{
    static StsObjectPanel objectPanel = null;
	static final long serialVersionUID = 1L;

    transient ForegroundThread fThread = null;
    transient String executeCmd = null;
    static boolean debug = false;
    transient float[] xyz = null;
    transient String exportDirectory = null;
    transient StsProgressBarDialog progressBarDialog = null;
    private transient StsMenuItem launchBtn = new StsMenuItem();
	private transient StsMenuItem exportBtn = new StsMenuItem();

    protected String commandString = null;
    protected String originalFilename = null;
    protected String filename = null;
    protected StsColor stsColor = new StsColor(StsColor.GREEN);
    protected double x = 0.0f;
    protected double y = 0.0f;
    protected double z = 0.0f;
    protected byte osDefinedOn = WINDOWS;
    protected StsWell well = null;

    static byte UNIX = 0;
    static byte WINDOWS = 1;
    static byte[] osTypes = { UNIX, WINDOWS };
    static String[] osTypesString = { "Unix", "Windows" };

    // Documents
    public static byte MSWORD = 0;
    public static byte MSXLS = 1;
    public static byte PDF = 2;
    public static byte TEXT = 3;
    public static byte[] doc_types = { MSWORD, MSXLS, PDF, TEXT };
    public static String[] doc_extensions = {"doc", "xls", "pdf", "txt"};
    public static String[] doc_msCommands = { "winword.exe", "excel.exe", "acrord32.exe", "notepad.exe" };
    public static String[] doc_unixCommands = { "", "", "", "" };

    // Images
    public static byte CGM = 0;
    public static byte JPG = 1;
    public static byte BMP = 2;
    public static byte GIF = 3;
    public static byte TIF = 4;
    public static byte[] img_types = { CGM, JPG, BMP, GIF, TIF };
    public static String[] img_extensions = {"cgm", "jpg", "bmp", "gif", "tif"};
    public static String[] img_msCommands = { "", "mspaint.exe", "mspaint.exe", "mspaint.exe", "mspaint.exe" };
    public static String[] img_unixCommands = { "", "", "", "" };

    // MultiMedia
    public static byte MPG = 0;
    public static byte AVI = 1;
    public static byte MP3 = 2;
    public static String[] mm_extensions = {"mpg", "avi", "mp3"};
    public static byte[] mm_types = { MPG, AVI, MP3 };
    public static String[] mm_msCommands = { "wmplayer.exe", "wmplayer.exe", "wmplayer.exe" };
    public static String[] mm_unixCommands = { "", "", "", "", "" };

    // Catch All Others
    public static byte OTHER = 0;
    public static String[] oth_extensions = {"*"};
    public static byte[] oth_types = { OTHER };
    public static String[] oth_msCommands = { "" };
    public static String[] oth_unixCommands = { "" };

    public static byte types[][] = {doc_types, img_types, mm_types, oth_types};
    public static String extensions[][] = {doc_extensions, img_extensions, mm_extensions, oth_extensions };
    public static String msCommands[][] = {doc_msCommands, img_msCommands, mm_msCommands, oth_msCommands };
    public static String unixCommands[][] = {doc_unixCommands, img_unixCommands, mm_unixCommands, oth_unixCommands };

    transient ImageIcon icon = null;

    static StsDateFieldBean bornField = new StsDateFieldBean(StsAncillaryData.class, "bornDate", "Born Date:");
    static StsDateFieldBean deathField = new StsDateFieldBean(StsAncillaryData.class, "deathDate", "Death Date:");

    static public StsFieldBean[] displayFields = null;

    public StsAncillaryData()
    {
    }

    public StsAncillaryData(String filename, String cmd, StsColor color, double x, double y, double z, boolean persist)
    {
        super(persist);

        try
        {
            this.originalFilename = filename;
            this.filename = currentModel.getProject().getDataFullDirString() +
                                       currentModel.getProject().getArchiveDirString() +
                                       StsFile.getFilenameFromPathname(filename);
            this.commandString = cmd;
            this.stsColor = color;
            this.x = x;
            this.y = y;
            this.z = z;

            icon = StsIcon.createIcon("dir16x32.gif");

            isVisible = true;

            setName(StsFile.getFilenameFromPathname(filename));
        }
        catch(Exception e) {}
    }

    public StsAncillaryData(String filename, StsColor color, boolean persist, byte type)
    {
        this(filename, getDefaultCommandFromFilename(filename, type), color, currentModel.getProject().getXOrigin(), currentModel.getProject().getYOrigin(),
             (double)currentModel.getProject().getZorTMin(), persist);
        setType(type);
    }

    public boolean initialize(StsModel model)
    {
        try
        {
            super.isVisible = isVisible;
            icon = StsIcon.createIcon("dir16x32.gif");
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsAncillaryData.classInitialize() failed.", e, StsException.WARNING);
            return false;
        }
    }
    public void display(StsGLPanel3d glPanel3d)
    {
        display(glPanel3d, false);
    }
    public void display(StsGLPanel3d glPanel3d, boolean displayNames)
    {
        if(!isVisible)
            return;

        try
        {
            float[] xy = currentModel.getProject().getRotatedRelativeXYFromUnrotatedAbsoluteXY(getXLoc(), getYLoc());
            xyz = new float[] { xy[0], xy[1], (float)getZLoc() };

            // Change this to an icon instead of a square...
            StsGLDraw.drawPoint(xyz, StsColor.BLACK, glPanel3d, 6, 1.0);
            StsGLDraw.drawPoint(xyz, stsColor, glPanel3d, 5, 2.0);

            if(displayNames)
            {
                stsColor.setGLColor(glPanel3d.getGL());
                GLBitmapFont font = GLHelvetica12BitmapFont.getInstance(glPanel3d.getGL());
                StsGLDraw.fontOutput(glPanel3d.getGL(), xyz, getName(), font);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsAncillaryData.display() failed.", e, StsException.WARNING);
        }
    }

    public void pick(GL gl, StsGLPanel glPanel)
    {
        StsGLDraw.drawPoint(xyz, gl, 6);
    }

    public void mouseSelectedEdit(StsMouse mouse)
    {
        reportMessage();
    }

    public void showPopupMenu(StsGLPanel glPanel, StsMouse mouse) { }
    public StsFieldBean[] getDisplayFields()
    {
        if(displayFields != null) return displayFields;
        {
            StsColorComboBoxFieldBean colorComboBoxBean =  new StsColorComboBoxFieldBean(StsAncillaryData.class, "stsColor", "Color:");
            StsSpectrum spectrum = currentModel.getSpectrum("Basic");
            colorComboBoxBean.setListItems(spectrum);
            displayFields = new StsFieldBean[]
            {
                new StsBooleanFieldBean(StsAncillaryData.class, "isVisible", "Enable:"),
                new StsStringFieldBean(StsAncillaryData.class, "name", true, "Name:"),
					colorComboBoxBean,
                new StsStringFieldBean(StsAncillaryData.class, "commandString", true, "Command:"),
                new StsDoubleFieldBean(StsAncillaryData.class, "xLoc", true, "X:"),
                new StsDoubleFieldBean(StsAncillaryData.class, "yLoc", true, "Y:"),
                new StsDoubleFieldBean(StsAncillaryData.class, "zLoc", true, "Z:"),
                new StsStringFieldBean(StsAncillaryData.class, "osAsString", false, "Defined on OS:")
           };
       }
       return displayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        return null;
    }

    public Object[] getChildren()
    {
        return new Object[0];
    }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }

    public boolean anyDependencies()
    {
        return false;
    }

    public void treeObjectSelected()
    {
        ;
    }

    public void setIsVisible(boolean b)
    {
        if (b == isVisible)
        {
            return;
        }
        isVisible = b;
        currentModel.win3dDisplayAll();
    }
    public boolean getIsVisible() { return isVisible; }

    public double getXLoc() { return x; }
    public double getYLoc() { return y; }
    public double getZLoc() { return z; }
    public StsPoint getLocation() { return new StsPoint(x,y,z); }
    public String getCommandString() { return commandString; }

    public void setXLoc(double value)
    {
        x = value;
        dbFieldChanged("x", x);
        currentModel.win3dDisplayAll();
    }
    public void setYLoc(double value)
    {
        y = value;
        dbFieldChanged("y", y);
        currentModel.win3dDisplayAll();
    }
    public void setZLoc(double value)
    {
        z = value;
        dbFieldChanged("z", z);
        currentModel.win3dDisplayAll();
    }
    public void setLocation(StsPoint point)
    {
        x = point.getX();
        y = point.getY();
        z = point.getZorT();
        currentModel.win3dDisplayAll();
    }
    public void setCommandString(String cmd) { commandString = cmd; }
    public void setStsColor(StsColor color)
    {
        if(this.stsColor == color) return;
        stsColor = color;
        currentModel.addTransactionCmd("ancillaryData color change", new StsChangeCmd(this, stsColor, "stsColor", false));
        currentModel.win3dDisplayAll();
    }

    static public byte getTypeFromFilename(String filename)
    {
        for(int j=0; j<extensions.length; j++)
        {
            for (int i = 0; i < extensions[j].length; i++)
            {
                if (filename.toLowerCase().endsWith(extensions[j][i]))
                {
                    return (byte)j;
                }
            }
        }
        return OTHER;
    }

    public String getOsAsString()
    {
        return osTypesString[osDefinedOn];
    }
    public void setOsAsString(String osName)
    {
        if(osName.indexOf("Windows") >= 0)
            osDefinedOn = WINDOWS;
        else
            osDefinedOn = UNIX;
    }
    public String getFilename() { return filename; }
    public String getOriginalFilename() { return originalFilename; }
    static public String[] getOSTypes() { return osTypesString; }
    public void updatePanel()
    {
        ;
    }

    public void reportMessage()
    {
        try
        {
            StsMessageFiles.infoMessage("Filename= " + StsFile.getFilenameFromPathname(filename));
        }
        catch(Exception e)
        {
            ;
        }
    }

    public StsWell getWell()
    {
        return well;
	}

    public void showPopup(MouseEvent e, StsGLPanel glPanel)
    {
        launchBtn.setMenuActionListener("Launch...", this, "launch", null);
		exportBtn.setMenuActionListener("Export...", this, "export", null);
        JPopupMenu tp = new JPopupMenu("Ancillary Data Popup");
        glPanel.add(tp);
        tp.add(launchBtn);
        tp.add(exportBtn);
        tp.show(currentModel.getGlPanel3d(), e.getX(), e.getY());
    }

    public void setWell(StsWell well)
    {
        this.well = well;
    }

    public boolean launch()
    {
        executeCmd = commandString + " " + filename;
        if(debug) System.out.println("Executing: " + executeCmd);
        if(fThread != null)
            fThread.killProcess();
        fThread = new ForegroundThread();
        fThread.start();
        return true;
    }

    class ForegroundThread extends Thread implements Runnable
   {
       /**
        * constructor. Automatically instantiates the object
        */
       Process prcs = null;
       ForegroundThread() { super(); }

       /**
        * Run the program
        */
       public void run()
       {
           try
           {
               Runtime rt = Runtime.getRuntime();
               prcs = rt.exec(executeCmd);
           }
           catch(Exception g)
           {
               new StsMessage(currentModel.win3d, StsMessage.ERROR, "Unable to launch: " + executeCmd +
                              "\n\n Verify the executable (" + commandString + ") is in your path.");
           }
           return;
       }

       public void killProcess()
       {
           if(prcs != null)
               prcs.destroy();
       }
   }

   public boolean defineLocation(StsGLPanel3d glPanel3d, StsMouse mouse)
   {
       if( mouse.getCurrentButton() != StsMouse.LEFT )
           return true;

       int buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
       if(buttonState == StsMouse.RELEASED)
       {
           StsProject project = currentModel.getProject();
           StsCursor3d cursor3d = currentModel.win3d.getCursor3d();
           if(cursor3d == null)
           {
               new StsMessage(currentModel.win3d, StsMessage.WARNING, "Must select point on cursor plane.");
               return false;
           }
           if(cursor3d.getCursorPoint(glPanel3d, mouse) == null)
               return false;

           StsCursorPoint cursorPoint = cursor3d.getNearestPointInCursorPlane(glPanel3d, mouse);
           StsPoint point = cursorPoint.point;
           if(point == null)
               return false;
           else
           {
               double[] xy = project.getAbsoluteXYCoordinates(point.getX(),point.getY());
               setXLoc(xy[0]);
               setYLoc(xy[1]);
               setZLoc(point.getZorT());

               currentModel.win3dDisplayAll();
           }
       }
       return true;
    }

    public boolean canExport() { return true; }
    public boolean canLaunch() { return true; }
    public boolean export()
    {
        try
        {
            StsDirectorySelectionDialog exportDialog = new StsDirectorySelectionDialog((Frame)currentModel.win3d,
                "Export Ancillary Data", StsFile.getFilenameFromPathname(filename), true);
            exportDialog.setVisible(true);
            exportDirectory = exportDialog.getCurrentDirectory();
            if(!exportDialog.getDoProcess())
            {
                return false;
            }
            progressBarDialog = StsProgressBarDialog.constructor(currentModel.win3d, "Export Ancillary Data", false);
            progressBarDialog.setLabelText("Exporting " + getName());
        }
        catch(Exception e)
        {
            new StsMessage(currentModel.win3d, StsMessage.ERROR,
                           "Problem launching directory selection screen for export");
        }

        Runnable runExport = new Runnable()
        {
            public void run()
            {
                exportData();
            }
        };

        Thread exportThread = new Thread(runExport);

        exportThread.start();
        return true;
    }

    private void exportData()
    {
        String outputName = null;
        try
        {
            outputName = exportDirectory + File.separator +
                StsFile.getFilenameFromPathname(filename);
            StsFile.copy(filename, outputName);
        }
        catch(Exception e)
        {
            new StsMessage(currentModel.win3d, StsMessage.ERROR, "Error exporting " + outputName);
        }
	}

    public void setBornDate(String born)
    {
        if(!StsDateFieldBean.validateDateInput(born))
        {
            bornField.setValue(StsDateFieldBean.convertToString(bornDate));
            return;
        }
        super.setBornDate(born);
    }

    public void setDeathDate(String death)
    {
        if(!StsDateFieldBean.validateDateInput(death))
        {
            deathField.setValue(StsDateFieldBean.convertToString(deathDate));
            return;
        }
        super.setDeathDate(death);
    }

	static public String getDefaultCommandFromFilename(String filename, byte type)
	{
		for(int i=0; i<extensions[type].length; i++)
		{
			if(filename.toLowerCase().endsWith(extensions[type][i]))
			{
				String osName = System.getProperty("os.name").toString();
				if(osName.indexOf("Windows") >= 0)
					return msCommands[type][i];
				else
					return unixCommands[type][i];
			}
		}
		return null;
	}

}
