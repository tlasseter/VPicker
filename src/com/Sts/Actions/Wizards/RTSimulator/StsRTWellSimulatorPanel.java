package com.Sts.Actions.Wizards.RTSimulator;

import com.Sts.Actions.Import.StsKeywordIO;
import com.Sts.Actions.Import.StsWellKeywordIO;
import com.Sts.Actions.Wizards.SensorLoad.StsSensorFile;
import com.Sts.DBTypes.StsMonitor;
import com.Sts.IO.StsAbstractFile;
import com.Sts.IO.StsAsciiFile;
import com.Sts.IO.StsFile;
import com.Sts.Interfaces.StsProgressRunnable;
import com.Sts.Types.StsWellFile;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Icons.StsIcon;
import com.Sts.UI.Progress.StsProgressPanel;
import com.Sts.UI.StsMessage;
import com.Sts.UI.StsYesNoDialog;
import com.Sts.Utilities.DateTime.CalendarParser;
import com.Sts.Utilities.StsException;
import com.Sts.Utilities.StsToolkit;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.StringTokenizer;

/*
 * ProgressBarDemo.java is a 1.4 application that requires these files:
 *   LongTask.java
 *   SwingWorker.java
 */

public class StsRTWellSimulatorPanel extends StsJPanel implements ActionListener {
    private String currentDirectory = "C:/Data/Demonstrations/ReservoirMonitoring/Realtime/Wells";
    private String currentDevFile = "C:/Data/Demonstrations/ReservoirMonitoring/Wells/well-dev.txt.1G-03RT";
    private String currentLogFile = null;
    private JFileChooser chooseDirectory = null;
    private JFileChooser chooseFile = null;
    private JFrame frame = null;
    private int rop = 25000; //feet per hour
    private String outputTypeString = StsMonitor.sourceStrings[StsMonitor.DIRECTORY];
    private boolean simulateFile = true;
    private StsIntFieldBean delayBean;

    protected OutputStream os = null;
    protected BufferedOutputStream bos = null;
    protected ByteArrayOutputStream baos = null;
    protected DataOutputStream dos = null;

    protected InputStream is = null;
    protected BufferedInputStream bis = null;
    protected DataInputStream dis = null;

    private processWellFiles devTask;
    private processWellFiles logTask;
    
    StsGroupBox inputBox = new StsGroupBox("Define Input (DD-MM-YY Ordered)");
    StsGroupBox internalBox = new StsGroupBox();
    StsGroupBox outputBox = new StsGroupBox("Define Output");
    StsGroupBox internalBox2 = new StsGroupBox();
    JLabel title = new JLabel();
    StsComboBoxFieldBean outputTypeBean = null;
    JTextField realtimeDirTxt = new JTextField();
    JTextField pathFileTxt = new JTextField();
    JTextField logFileTxt = new JTextField();
    StsDateFieldBean dateBean = new StsDateFieldBean();
    JButton directoryBrowseButton = new JButton();
    JButton fileDevBrowseButton = new JButton();
    JButton fileLogBrowseButton = new JButton();
    JButton acceptBtn = new JButton();
    JButton cancelBtn = new JButton();
    JButton exitBtn = new JButton();

    SimpleDateFormat format = new SimpleDateFormat("dd/MM/yy kk:mm:ss.SSS");
    SimpleDateFormat outFormat = new SimpleDateFormat("dd/MM/yy kk:mm:ss.SSS");
    String startDate = new String("01/01/00 00:00:00.000");

	public StsProgressPanel progressPanel = StsProgressPanel.constructorWithCancelButton(5, 50);

    public StsRTWellSimulatorPanel()
    {
        dateBean = new StsDateFieldBean(this, "startDate", "01/01/01 00:00:00.000", true, "Start: ");

        outputTypeBean = new StsComboBoxFieldBean(this, "outputTypeString", "Source Type: ", StsMonitor.sourceStrings);
        outputTypeBean.setEditable(false);

        devTask = new processWellFiles(progressPanel);
        logTask = new processWellFiles(progressPanel);

        delayBean = new StsIntFieldBean(this, "delay", 100, 100000, "Rate of Penetration (ft/hr):", true);
        delayBean.setValue(rop);
        delayBean.setToolTipText("Increase the speed of data flow by increasing the rate of panetration.");

        fileDevBrowseButton.setText("Deviation File");
        fileDevBrowseButton.setBorder(BorderFactory.createRaisedBevelBorder());
        pathFileTxt.setText("<supply well deviation file>");
        pathFileTxt.setBackground(Color.WHITE);
        pathFileTxt.setToolTipText("Press selection button or type input file name followed by return/enter key.");
        fileLogBrowseButton.setText("Log File");
        fileLogBrowseButton.setBorder(BorderFactory.createRaisedBevelBorder());
        logFileTxt.setText("<supply well log file>");
        logFileTxt.setBackground(Color.WHITE);
        logFileTxt.setToolTipText("Press selection button or type input file name followed by return/enter key.");

        directoryBrowseButton.setIcon(StsIcon.createIcon("dir16x32v2.gif"));
        directoryBrowseButton.setBorder(BorderFactory.createRaisedBevelBorder());
        realtimeDirTxt.setText("<supply output directory>");
        realtimeDirTxt.setBackground(Color.WHITE);
        realtimeDirTxt.setToolTipText("Press selection button or type output dirNo or file name followed by return/enter key.");


        acceptBtn.setText("Accept");
        cancelBtn.setText("Cancel");
        exitBtn.setText("Exit");

        internalBox.gbc.weightx = 0.0;
        internalBox.gbc.gridwidth = 1;
        internalBox.addToRow(fileDevBrowseButton);
        internalBox.gbc.weightx = 1.0;
        internalBox.gbc.gridwidth = 4;
        internalBox.gbc.fill = inputBox.gbc.HORIZONTAL;
        internalBox.addEndRow(pathFileTxt);
        internalBox.gbc.weightx = 0.0;
        internalBox.gbc.gridwidth = 1;
        internalBox.addToRow(fileLogBrowseButton);
        internalBox.gbc.weightx = 1.0;
        internalBox.gbc.gridwidth = 4;
        internalBox.gbc.fill = inputBox.gbc.HORIZONTAL;
        internalBox.addEndRow(logFileTxt);

        inputBox.gbc.fill = inputBox.gbc.HORIZONTAL;
        inputBox.gbc.gridwidth = 5;
        inputBox.addEndRow(internalBox);
        inputBox.addToRow(dateBean);

        outputBox.gbc.fill = inputBox.gbc.HORIZONTAL;
        outputBox.addEndRow(outputTypeBean);

        internalBox2.addToRow(directoryBrowseButton);
        internalBox2.gbc.fill = internalBox2.gbc.HORIZONTAL;
        internalBox2.addEndRow(realtimeDirTxt);

        outputBox.addEndRow(internalBox2);
        outputBox.addEndRow(delayBean);

        gbc.fill = gbc.HORIZONTAL;
        gbc.anchor = NORTH;
        gbc.weighty = 0.0;
        addEndRow(inputBox);
        addEndRow(outputBox);
        gbc.weighty = 1.0;
        gbc.fill = gbc.BOTH;        
        addEndRow(progressPanel);
        gbc.fill = gbc.HORIZONTAL;
        gbc.weighty = 0.0;        
        addEndRow(acceptBtn);

        directoryBrowseButton.addActionListener(this);
        fileDevBrowseButton.addActionListener(this);
        fileLogBrowseButton.addActionListener(this);
        acceptBtn.addActionListener(this);
        realtimeDirTxt.addActionListener(this);

        setCurrentDirectory(currentDirectory);
        setCurrentDevFile(currentDevFile);
        setCurrentLogFile(currentLogFile);
    }

    /**
     * Called when the user presses the start button.
     */
    public void actionPerformed(ActionEvent evt) {
        Object source = evt.getSource();
        if (source == fileDevBrowseButton)
        {
            chooseFile = new JFileChooser(currentDevFile);
            chooseFile.setCurrentDirectory(new File(currentDirectory));
            chooseFile.setFileSelectionMode(JFileChooser.FILES_ONLY);
            chooseFile.setDialogTitle("Select the input deviation file.");
            chooseFile.setApproveButtonText("Open File");
            while(true)
            {
                chooseFile.showOpenDialog(null);
                if(chooseFile.getSelectedFile() == null)
                    break;
                File newFile = chooseFile.getSelectedFile();
                if (newFile.isFile())
                {
                    setCurrentDevFile(newFile.getAbsolutePath());
                    break;
                }
                else
                {
                    if(!StsYesNoDialog.questionValue(this,"Must select the file that\n contains the sensor data.\n\n Continue?"))
                        break;
                }
            }
        }
        else if (source == fileLogBrowseButton)
        {
            chooseFile = new JFileChooser(currentLogFile);
            chooseFile.setCurrentDirectory(new File(currentDirectory));
            chooseFile.setFileSelectionMode(JFileChooser.FILES_ONLY);
            chooseFile.setDialogTitle("Select the input log file.");
            chooseFile.setApproveButtonText("Open File");
            while(true)
            {
                chooseFile.showOpenDialog(null);
                if(chooseFile.getSelectedFile() == null)
                {
                    currentLogFile = null;
                    break;
                }
                File newFile = chooseFile.getSelectedFile();
                if (newFile.isFile())
                {
                    setCurrentLogFile(newFile.getAbsolutePath());
                    break;
                }
                else
                {
                    if(!StsYesNoDialog.questionValue(this,"Must select the file that\n contains the log data.\n\n Continue?"))
                        break;
                }
            }
        }
        else if (source == directoryBrowseButton)
        {
            chooseDirectory = new JFileChooser(currentDirectory);
            chooseDirectory.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            chooseDirectory.setDialogTitle("Select or Enter Desired File or Directory and Press Open");
            chooseDirectory.setApproveButtonText("Open");
            chooseDirectory.showOpenDialog(null);
            File newDirFile = chooseDirectory.getSelectedFile();
            if(newDirFile != null)
            	setCurrentDirectory(newDirFile.getAbsolutePath());
        }
        else if(source == realtimeDirTxt)  // Need carriage return to accept
        {
        	setCurrentDirectory(realtimeDirTxt.getText());
        }
        else if(source == acceptBtn)
        {
        	if(devTask.done == true)
        	{
        		if(!validateUserInput())
        			return;
        		StsProgressRunnable progressRunnable = devTask.getProcessDevFileRunnable();
                Thread devthread = new Thread(progressRunnable);
                devthread.start();
        		//StsToolkit.runRunnable(progressRunnable);
                StsProgressRunnable progressLogRunnable = logTask.getProcessLogFileRunnable();
                Thread logthread = new Thread(progressLogRunnable);
                logthread.start();
        		//StsToolkit.runRunnable(progressRunnable);
        	}
        	else
        	{
        		if(StsYesNoDialog.questionValue(frame, "Do you want to re-start?"))
        		{
        			devTask.progressPanel.cancel();
        			while(!devTask.done)
        				try {Thread.sleep(rop);} catch(Exception e) { }

            		StsProgressRunnable progressRunnable = devTask.getProcessDevFileRunnable();
            		StsToolkit.runRunnable(progressRunnable);
                    progressRunnable = logTask.getProcessLogFileRunnable();
        		    StsToolkit.runRunnable(progressRunnable);
        		}
        	}
        }
    }

    public void setOutputTypeString(String typeString)
    {
    	if(typeString.equals(StsMonitor.sourceStrings[StsMonitor.DIRECT]) || typeString.equals(StsMonitor.sourceStrings[StsMonitor.DB]))
    	{
    		new StsMessage(this.frame, StsMessage.WARNING, typeString + " is currently not supported, pick a different type.");
    		setOutputTypeString(StsMonitor.sourceStrings[StsMonitor.DIRECTORY]);
    	}
    	outputTypeString = typeString;
    	if(outputTypeString.equalsIgnoreCase(StsMonitor.sourceStrings[StsMonitor.DIRECTORY]))
    		directoryBrowseButton.setIcon(StsIcon.createIcon("dir16x32v2.gif"));
    	else
    		directoryBrowseButton.setIcon(StsIcon.createIcon("file16x32v2.gif"));
    	this.repaint();
    }
    public String getOutputTypeString() { return outputTypeString; }
    public boolean validateUserInput()
    {
    	File fileDir = null;
    	try
    	{
            if(!currentDevFile.equalsIgnoreCase("None"))
            {
    	        fileDir = new File(currentDevFile);
    		    if(!fileDir.isFile())
    		    {
    			    new StsMessage(this.frame, StsMessage.ERROR, "Selected file is not valid file.");
    			    return false;
    		    }
            }

        	fileDir = new File(currentDirectory);
        	if(fileDir.isFile())
                currentDirectory = fileDir.getAbsolutePath();
            return true;
    	}
    	catch(Exception e)
    	{
    		new StsMessage(this.frame, StsMessage.ERROR, "Selected input file and/or output directory are not valid.");
    		return false;
    	}
    }
    public int getDelay() { return rop; }
    public void setDelay(int delay) { this.rop = delay; }
    public boolean getSimulateFile() { return simulateFile; }
    public void setSimulateFile(boolean time)
    {
    	this.simulateFile = time;
    }
    private void setCurrentDirectory(String directory)
    {
        currentDirectory = directory;
        realtimeDirTxt.setText(currentDirectory);
        validateUserInput();
    }
    public String getStartDate() { return startDate; }
    public void setStartDate(String start) { startDate = start; }
    private void setCurrentDevFile(String filename)
    {
        currentDevFile = filename;
        pathFileTxt.setText(currentDevFile);
    }
    private void setCurrentLogFile(String filename)
    {
        currentLogFile = filename;
        logFileTxt.setText(currentLogFile);
    }

private class processWellFiles {
    private boolean done = true;
    private boolean canceled = false;
    private StsProgressPanel progressPanel;

    public processWellFiles(StsProgressPanel progressPanel) {
    	this.progressPanel = progressPanel;
    }

    private StsProgressRunnable getProcessDevFileRunnable()
    {
        StsProgressRunnable progressRunnable = new StsProgressRunnable()
        {
            public void cancel()
            {
                stop();
            }
            public void run()
            {
            	done = false;
            	progressPanel.progressBar.canceled = false;
                try
                {
                    StsKeywordIO.parseAsciiFilename(StsFile.getFilenameFromPathname(currentDevFile));
                    String group = StsKeywordIO.group;
                    String stubname = StsKeywordIO.getFileStemName(currentDevFile);
                    readFile(currentDevFile, stubname, group);
                }
                catch(Exception ex)
                {
                    StsException.printStackTrace(ex);
                }
            }
        };
        return progressRunnable;
    }

    private StsProgressRunnable getProcessLogFileRunnable()
    {
        StsProgressRunnable progressRunnable = new StsProgressRunnable()
        {
            public void cancel()
            {
                stop();
            }
            public void run()
            {
                if(currentLogFile == null)
                    return;
            	done = false;
            	progressPanel.progressBar.canceled = false;
                try
                {
                    StsKeywordIO.parseAsciiFilename(StsFile.getFilenameFromPathname(currentLogFile));
                    String group = StsKeywordIO.group;
                    String stubname = StsKeywordIO.getFileStemName(currentLogFile);
                    readFile(currentLogFile, stubname, group);
                }
                catch(Exception ex)
                {
                    StsException.printStackTrace(ex);
                }
            }
        };
        return progressRunnable;
    }

    public void stop()
    {
        canceled = true;
        progressPanel.progressBar.canceled = true;
        done = true;
    }

    /**
     * Called from ProgressBarDemo to find out if the task has completed.
     */
    public boolean isDone()
    {
        return done;
    }

    public void readFile(String inputFile, String stubname, String group)
    {
    	StsFile file = null;
        StsAsciiFile aFile = null;

        String line = null;
        file = StsFile.constructor(inputFile);

        aFile = new StsAsciiFile(file);
        aFile.openRead();
        //System.out.println("Input file= " + inputFile);
        try
        {
            long previousTime = ((Date)format.parse(startDate)).getTime() - 1;
            long currentTime = previousTime;

            int current = 1;
            StsAbstractFile oFile = null;
            StsAsciiFile outFile = null;

            String oFilename = currentDirectory;
            //
            // If creating file, add header once and then loop through lines.
            //
            StsWellFile wellFile = new StsWellFile(file);
            wellFile.analyzeFile();
            int numLines = wellFile.getNumLinesInFile();
            if(group.equalsIgnoreCase("well-dev"))
            {
                oFilename = oFilename + "\\well-dev.txt." + stubname;
                outFile = createDevOutputFile(oFilename, stubname, wellFile);
                progressPanel.setMaximum(numLines);
            }
            else
            {
                oFilename = oFilename + "\\well-logs.txt." + stubname;
                outFile = createLogOutputFile(oFilename, stubname, wellFile);
            }
            int depthIdx = wellFile.getColLocation(wellFile.Z); // Might be depth, would prefer mdepth
            for(int j=0; j< wellFile.curveNames.length; j++)
            {
                if(wellFile.curveNames[j].equalsIgnoreCase("MDEPTH"))
                    depthIdx = j;
            }
            //System.out.println("Output file= " + oFilename);
            float previousDepth = 0.0f;
            wellFile.analyzeFile();
            while(!done)
            {
                // Read or create event
                line = wellFile.readLine();
                if(line == null)
                   	break;

                // Compute the change in depth
                String token = null;
                float depth = 0.0f;
                StringTokenizer stok = new StringTokenizer(line," \t");
                String[] tokens = new String[stok.countTokens()];
                for(int i=0; i<tokens.length; i++)
                {
                    token = stok.nextToken();
                    if(depthIdx == i)
                        depth = Float.valueOf(token).floatValue();
                }
                float depthDiff = depth - previousDepth;
                previousDepth = depth;

                // Wait the requested delay.
                double ropMs = rop/3600000.0f; // Feet per millisecond
                currentTime = previousTime + Math.abs((long)(depthDiff/ropMs));
                long waitTime = currentTime - previousTime;
                if(waitTime < 0)
                   waitTime = 10000;
                //System.out.println("Waiting " + waitTime + " ms on line: " + line);
                Thread.sleep(waitTime);
                previousTime = currentTime;
                // Output line
                // --- If simulate real-time, replace time in file with current time.
                long time = System.currentTimeMillis();
                if(group.equalsIgnoreCase("well-dev"))
                {
                    String timeDate = outFormat.format(new Date(time));
                    line = timeDate + " " + line;
                }
                //System.out.println("Outputting " + line + " to file: " + outFile.getFilename());
                boolean writeSuccessful = false;
                while(!writeSuccessful)
                {
                    try
                    {
                        outFile.openWriteAppend();
                        outFile.writeLine(line);
                        if(group.equalsIgnoreCase("well-dev"))
               	            progressPanel.setValue(current-1);
               	        if(StsMonitor.getTypeFromString(outputTypeString) == StsMonitor.DIRECTORY)
                        {
                   	        progressPanel.appendLine("Output file for time " + format.format(new Date(time)) + " at depth " + depth + " to " + outFile.getFilename());
                        }
                        else
                   	        progressPanel.appendLine("Output line #" + current + " at depth " + depth + " at time " + format.format(new Date(time)) + " to file: " + outFile.getFilename());
                        current++;

               	        if(progressPanel.isCanceled())
                    	{
               		        progressPanel.appendLine("Canceling processing... ");
               		        done = true;
               	        }
                        outFile.close();
                        writeSuccessful = true;
                    }
                    catch(Exception ex)
                    {
                        continue;
                    }
                }

            }
            if(!canceled)
            {
               	progressPanel.finished();
               	progressPanel.setDescriptionAndLevel("Processing complete", progressPanel.progressBar.INFO);
            }
            else
               	progressPanel.setDescriptionAndLevel("Canceled by user", progressPanel.progressBar.WARNING);
            }
            catch (Exception e) {
            	 progressPanel.setDescriptionAndLevel("Error processing file. Check permissions.", progressPanel.progressBar.ERROR);
                 return;
             }
             finally
             {
                 if(!simulateFile)
                 {
                     aFile.close();
                 }
                 done = true;
             }
          }
       }

    private StsAsciiFile createDevOutputFile(String filename, String stubname, StsWellFile wFile)
    {
        StsAsciiFile outfile = null;
        StsAbstractFile oFile = null;

        try
        {
            oFile = StsFile.constructor(filename);
            outfile = new StsAsciiFile(oFile);

            // Truncate file and start on line one.
            outfile.openWrite();
            progressPanel.appendLine("Output header to file: " + outfile.getFilename());
            outfile.writeLine("WELLNAME");
            outfile.writeLine(stubname);
            outfile.writeLine("ORIGIN XY");
            outfile.writeLine(wFile.xOrigin + " " + wFile.yOrigin);
            outfile.writeLine("CURVE");
            outfile.writeLine("DATE");
            outfile.writeLine("CTIME");
            for(int i=0; i<wFile.curveNames.length; i++)
                outfile.writeLine(wFile.curveNames[i]);
            outfile.writeLine("VALUE");
            outfile.close();
            return outfile;
        }
        catch(Exception ex)
        {
            return null;
        }
    }


    private StsAsciiFile createLogOutputFile(String filename, String stubname, StsWellFile wFile)
    {
        StsAsciiFile outfile = null;
        StsAbstractFile oFile = null;

        try
        {
            oFile = StsFile.constructor(filename);
            outfile = new StsAsciiFile(oFile);

            // Truncate file and start on line one.
            outfile.openWrite();
            progressPanel.appendLine("Output header to file: " + outfile.getFilename());
            outfile.writeLine("WELLNAME");
            outfile.writeLine(stubname);
            outfile.writeLine("CURVE");
            for(int i=0; i<wFile.curveNames.length; i++)
                outfile.writeLine(wFile.curveNames[i]);
            outfile.writeLine("VALUE");
            outfile.close();
            return outfile;
        }
        catch(Exception ex)
        {
            return null;
        }
    }
}