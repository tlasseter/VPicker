package com.Sts.Actions.Wizards.KMSTranslator;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;

/*
 * ProgressBarDemo.java is a 1.4 application that requires these files:
 *   LongTask.java
 *   SwingWorker.java
 */
public class StsKMSTranslatorDialog extends JPanel implements ActionListener {
    private String currentFile = "C:/Data/Projects/KMSTech";
    private JFileChooser chooseFile = null;    
    private JFrame frame = null;
    private boolean simulateFile = false;
    
    protected String ofilename = "aramco";
    protected OutputStream os = null;
    protected BufferedOutputStream bos = null;
    protected ByteArrayOutputStream baos = null;
    protected DataOutputStream dos = null;

    protected InputStream is = null;
    protected BufferedInputStream bis = null;
    protected DataInputStream dis = null;

    protected int numLines = 0;
    protected float xMin, xMax, xInc = -1.0f, yMin, yMax, yInc = -1.0f, zMin, zMax, zInc = -1.0f, aMin, aMax;
    protected int numInlines = -1, numXLines = -1, numSlices = -1;
    protected boolean analyze = false;
    protected float[] attribute;

    private processFile task;

    JPanel jPanel1 = new JPanel();
    JPanel jPanel2 = new JPanel();
    StsGroupBox inputBox = new StsGroupBox("Define Input");
    StsGroupBox fileBox = new StsGroupBox("File Statistics");
    StsFloatFieldBean xMinBean, xMaxBean, yMinBean, yMaxBean, zMinBean, zMaxBean, aMinBean, aMaxBean;
    StsIntFieldBean numInlinesBean, numXlinesBean, numSlicesBean;

    JPanel mainPanel = new JPanel();
    JLabel title = new JLabel();
    StsComboBoxFieldBean outputTypeBean = null;    
    JTextField sensorFileTxt = new JTextField();
    JButton fileBrowseButton = new JButton();
    JButton acceptBtn = new JButton();
    JButton cancelBtn = new JButton();
    JButton exitBtn = new JButton();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    GridBagLayout gridBagLayout3 = new GridBagLayout();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout4 = new GridBagLayout();
    
	public StsProgressPanel progressPanel = StsProgressPanel.constructorWithCancelButton(5, 50);
	
    public StsKMSTranslatorDialog()
    {
        super(new BorderLayout());

        xMinBean = new StsFloatFieldBean(this, "xMin", false, "X Minimum:");
        xMaxBean = new StsFloatFieldBean(this, "xMax", false, "X Maximum:");
        yMinBean = new StsFloatFieldBean(this, "yMin", false, "Y Minimum:");
        yMaxBean = new StsFloatFieldBean(this, "yMax", false, "Y Maximum:");
        zMinBean = new StsFloatFieldBean(this, "zMin", false, "Z Minimum:");
        zMaxBean = new StsFloatFieldBean(this, "zMax", false, "Z Maximum:");
        aMinBean = new StsFloatFieldBean(this, "aMin", false, "Resistivity Minimum:");
        aMaxBean = new StsFloatFieldBean(this, "aMax", false, "Resistivity Maximum:");
        numInlinesBean = new StsIntFieldBean(this, "numInlines", false, "Num Inlines:");
        numXlinesBean = new StsIntFieldBean(this, "numXLines", false, "Num XLines:");
        numSlicesBean = new StsIntFieldBean(this, "numSlices", false, "Num Slices:");

        task = new processFile(progressPanel);
        this.setLocation(200,200);
        jPanel1.setBorder(BorderFactory.createEtchedBorder());
        jPanel1.setLayout(gridBagLayout4);
        jPanel1.setBackground(SystemColor.menu);
        jPanel2.setBorder(BorderFactory.createEtchedBorder());
        jPanel2.setLayout(gridBagLayout2);
        jPanel2.setBackground(SystemColor.menu);
        title.setFont(new java.awt.Font("KMS Translator", 1, 14));
        title.setForeground(new Color(0, 55, 152));
        title.setOpaque(false);
        title.setHorizontalAlignment(SwingConstants.CENTER);
        title.setText("Convert KMS EM Data to S2S Format");
        sensorFileTxt.setText("<supply kms file>");
        sensorFileTxt.setBackground(Color.WHITE);   
        sensorFileTxt.setToolTipText("Press selection button or type input file name followed by return/enter key.");
        fileBrowseButton.setIcon(StsIcon.createIcon("file16x32v2.gif"));
        fileBrowseButton.setBorder(BorderFactory.createRaisedBevelBorder());
        acceptBtn.setText("Accept");
        cancelBtn.setText("Cancel");
        exitBtn.setText("Exit");
        acceptBtn.setBackground(SystemColor.menu);
        cancelBtn.setBackground(SystemColor.menu);
        exitBtn.setBackground(SystemColor.menu);

        mainPanel.setLayout(gridBagLayout1);
        mainPanel.setBackground(SystemColor.menu);
        mainPanel.setMinimumSize(new Dimension(250, 300));
        mainPanel.setPreferredSize(new Dimension(700, 500));

        inputBox.gbc.weightx = 0.0;   
        inputBox.gbc.gridwidth = 1;
        inputBox.addToRow(fileBrowseButton);
        inputBox.gbc.weightx = 1.0; 
        inputBox.gbc.gridwidth = 4;
        inputBox.gbc.fill = inputBox.gbc.HORIZONTAL;
        inputBox.addEndRow(sensorFileTxt);

        fileBox.gbc.fill = fileBox.gbc.HORIZONTAL;
        fileBox.addToRow(xMinBean);
        fileBox.addToRow(xMaxBean);
        fileBox.addEndRow(numXlinesBean);
        fileBox.addToRow(yMinBean);
        fileBox.addToRow(yMaxBean);
        fileBox.addEndRow(numInlinesBean);
        fileBox.addToRow(zMinBean);
        fileBox.addToRow(zMaxBean);
        fileBox.addEndRow(numSlicesBean);
        fileBox.addToRow(aMinBean);
        fileBox.addEndRow(aMaxBean);

        mainPanel.add(jPanel1,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 0, 0), 5, 6));

        jPanel1.add(title,  new GridBagConstraints(0, 0, 2, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 4, 0, 5), 181, 11));
        jPanel1.add(inputBox,  new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(2, 4, 0, 5), 0, 0));
        jPanel1.add(fileBox,  new GridBagConstraints(0, 2, 1, 1, 1.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(2, 4, 0, 5), 0, 0));
        jPanel1.add(progressPanel,   new GridBagConstraints(0, 5, 2, 1, 0.0, 1.0
                ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 16, 5, 5), 0, 0));
        
        mainPanel.add(jPanel2,  new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 5, 5, 0), 250, 0));
        jPanel2.add(acceptBtn, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 0, 0));

        add(mainPanel);

        fileBrowseButton.addActionListener(this);
        acceptBtn.addActionListener(this);
        sensorFileTxt.setText(currentFile);
    }

    /**
     * Called when the user presses the start button.
     */
    public void actionPerformed(ActionEvent evt) {
        Object source = evt.getSource();
        if (source == fileBrowseButton)
        {
            chooseFile = new JFileChooser(currentFile);
            chooseFile.setFileSelectionMode(JFileChooser.FILES_ONLY);
            chooseFile.setDialogTitle("Select the input sensor file.");
            chooseFile.setApproveButtonText("Open File");
            while(true)
            {
                chooseFile.showOpenDialog(null);
                if(chooseFile.getSelectedFile() == null)
                    break;
                File newFile = chooseFile.getSelectedFile();
                if (newFile.isFile())
                {
                    setCurrentFile(newFile.getAbsolutePath());
                    break;
                }
                else
                {
                    if(!StsYesNoDialog.questionValue(this,"Must select KMS EM file.\n\n Continue?"))
                        break;
                }
            }
        }
        else if(source == acceptBtn)
        {       	
        	if(task.done == true)
        	{
                analyze = false;
                progressPanel.initialize(numLines);
        		StsProgressRunnable progressRunnable = task.getProcessFileRunnable();
        		StsToolkit.runRunnable(progressRunnable);
        	}
        }
    }

    private void setCurrentFile(String filename)
    {
        analyze = true;
        currentFile = filename;
        sensorFileTxt.setText(currentFile);        
        if(task.done == true)
        {
        	StsProgressRunnable progressRunnable = task.getProcessFileRunnable();
        	StsToolkit.runRunnable(progressRunnable);
        }

    }

    /**
     * Create the GUI and show it.  For thread safety,
     * this method should be invoked from the
     * event-dispatching thread.
     */
    private static void createAndShowGUI() {
        JFrame frame = new JFrame("KMS File Translator");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        //Create and set up the content pane.
        JComponent newContentPane = new StsKMSTranslatorDialog();
        newContentPane.setOpaque(true); //content panes must be opaque
        frame.setContentPane(newContentPane);

        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }
    
    public static void main(String[] args)
    {
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { createAndShowGUI(); }});
    }
    
private class processFile {
    private boolean done = true;
    private boolean canceled = false;
    private StsProgressPanel progressPanel;
    
    public processFile(StsProgressPanel progressPanel) {
    	this.progressPanel = progressPanel;
    }
    
    private StsProgressRunnable getProcessFileRunnable()
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
                if(analyze == false)
                    readFile();
                else
                    analyzeFile();
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

    private void analyzeFile()
    {
    	int cnt = 0;
        float x, y, z, a;
        String line = null;
    	StsAbstractFile iFile = StsFile.constructor(currentFile);
        int fileSize = (int)iFile.getFile().length();
        progressPanel.initialize(fileSize);
        int estimatedLines = 0;
        StsAsciiFile aiFile = new StsAsciiFile(iFile);
        try
        {
            aiFile.openRead();
            int bytesRead = 0;
            progressPanel.setDescription("Analyzing File Ranges....Stand by");
            float firstX = -1.0f, firstY = -1.0f, firstZ = -1.0f;
            while((line = aiFile.readLine()) != null)
            {
                bytesRead += line.length();
                if(numLines == 0)
                {
                    estimatedLines = fileSize/bytesRead;
                    attribute = new float[estimatedLines+100];
                }
                progressPanel.setValue(bytesRead);
				line = line.trim();
                line = StsStringUtils.detabString(line);

                // Determine inline, crossline and slice ranges
			    StringTokenizer stok = new StringTokenizer(line);
 				x = Float.valueOf(stok.nextToken()).floatValue();
                if(xMin > x) xMin = x;
                if(xMax < x) xMax = x;
                if(numLines == 0)
                    firstX = x;
                else
                {
                    if((firstX != x) && (xInc == -1.0f))
                        xInc = Math.abs(Math.abs(firstX) - Math.abs(x));
                }
				y = Float.valueOf(stok.nextToken()).floatValue();
                if(yMin > y) yMin = y;
                if(yMax < y) yMax = y;
                if(numLines == 0)
                    firstY = y;
                else
                {
                    if((firstY != y) && (yInc == -1.0f))
                        yInc = Math.abs(Math.abs(firstY) - Math.abs(y));
                }
				z = -Float.valueOf(stok.nextToken()).floatValue();   // Convert to depth from elevation.
                if(zMin > z) zMin = z;
                if(zMax < z) zMax = z;
                if(numLines == 0)
                    firstZ = z;
                else
                {
                    if((firstZ != z) && (zInc == -1.0f))
                        zInc = Math.abs(Math.abs(firstZ) - Math.abs(z));
                }
				a = Float.valueOf(stok.nextToken()).floatValue();
                attribute[numLines] = a;
                if(aMin > a) aMin = a;
                if(aMax < a) aMax = a;

        	    numLines++;
                if((numLines%1000) == 0) progressPanel.setDescription("Processed " + numLines + " lines...");
                if((numLines%100000) == 0) progressPanel.appendLine("Processed " + numLines + " lines...");
            }
            attribute = (float[])StsMath.trimArray(attribute, numLines);
            // Compute the number of inlines, crosslines and slices.
            numSlices = (int)((zMax-zMin)/zInc) + 1;
            numXLines = (int)((xMax-xMin)/xInc) + 1;
            numInlines = (int)((yMax-yMin)/yInc) + 1;
            progressPanel.setValue(fileSize);
            progressPanel.setDescription("Analysis complete..." + numLines + " lines read.");
            xMinBean.setValue(xMin);
            xMaxBean.setValue(xMax);
            numInlinesBean.setValue(numInlines);
            yMinBean.setValue(yMin);
            yMaxBean.setValue(yMax);
            numXlinesBean.setValue(numXLines);
            zMinBean.setValue(zMin);
            zMaxBean.setValue(zMax);
            numSlicesBean.setValue(numSlices);
            aMinBean.setValue(aMin);
            aMaxBean.setValue(aMax);
        }
        catch(Exception e)
        {
        	System.out.println("Failed to count number of lines in file:" + currentFile + " Exception=" + e);
        }
        finally
        {
            done = true;
            aiFile.close();
        }
        return;
    }

    public void readFile()
    {
    	StsAbstractFile iFile = null;
        StsAsciiFile aiFile = null;
        String prefix = "seis3d.";
        iFile = StsFile.constructor(currentFile);
        aiFile = new StsAsciiFile(iFile);

        StsSegyVolume volume = new StsSegyVolume();
        volume.setOrigin(xMin, yMin);        
        volume.setAngle(0.0f);
        volume.setAngleSet(true);
        volume.setColNumInc(1.0f);
        volume.setColNumMax(numXLines);
        volume.setColNumMin(0);
        volume.setRowNumInc(1.0f);
        volume.setRowNumMax(numInlines);
        volume.setRowNumMin(0);
        volume.setZInc(zInc);
        volume.setZMax(zMax);
        volume.setZMin(zMin);
        volume.setXInc(xInc);
        volume.setXMax(xMax);
        volume.setXMin(xMin);
        volume.setYInc(yInc);
        volume.setYMax(yMax);
        volume.setYMin(yMin);
        volume.setZDomainString(StsParameters.TD_DEPTH_STRING);
        volume.setDataMax(aMax);
        volume.setDataMin(aMin);
        volume.setIsRegular(true);
        volume.setName("KMS");
        volume.setNCols(numXLines);
        volume.setNRows(numInlines);
        volume.setNSlices(numSlices);
        int dotLoc = currentFile.lastIndexOf(".");
        int slashLoc = currentFile.lastIndexOf("\\");
        volume.setRowCubeFilename(prefix + "bytes-inline." + currentFile.substring(slashLoc+1,dotLoc));
        volume.setColCubeFilename(prefix + "bytes-xline." + currentFile.substring(slashLoc+1,dotLoc));
        volume.setSliceCubeFilename(prefix + "bytes-trace." + currentFile.substring(slashLoc+1,dotLoc));
        try
        {
            StsAbstractFile oFile = null;
            StsBinaryFile outfile = null;
            float onePrct = (float)(aMax-aMin)*.01f;
            // Output Header File
            StsParameterFile.writeObjectFields(currentFile.substring(0, slashLoc+1) + prefix + "txt." + currentFile.substring(slashLoc+1,dotLoc)
                    ,volume, StsSeismicBoundingBox.class, StsMainObject.class);
            progressPanel.appendLine("Output header file: " + currentFile.substring(0, slashLoc+1) + prefix + "txt." + currentFile.substring(slashLoc+1,dotLoc));

            // Output Inline  -- Need to sort
            String fileName = currentFile.substring(0, slashLoc+1) + "seis3d.bytes-inline." + currentFile.substring(slashLoc+1,dotLoc);
            progressPanel.setMaximum(3);
            progressPanel.appendLine("Outputting inline file: " + fileName);
            FileOutputStream fileoutputstream = new FileOutputStream(fileName);
            DataOutputStream dataoutputstream = new DataOutputStream(fileoutputstream);
            byte[] bytes = null;
            float[] attributeLine = new float[numSlices];
            for(int i=0; i<numInlines; i++)
            {
                for(int n=0; n<numXLines; n++)
                {
                    int offset = i*numSlices + n*numInlines*numSlices;
                    System.arraycopy(attribute, offset, attributeLine, 0, numSlices);                                        
                    bytes = StsMath.floatsToUnsignedBytes254(attributeLine, aMin, aMax+onePrct);
                    dataoutputstream.write(bytes);
                }
            }
            progressPanel.appendLine("Completed Output of inline file: " + fileName);
            progressPanel.setValue(1);

            // Output Crossline  -- Should be okay
            bytes = null;
            fileName = currentFile.substring(0, slashLoc+1) + "seis3d.bytes-xline." + currentFile.substring(slashLoc+1,dotLoc);
            progressPanel.appendLine("Outputting crossline file: " + fileName);
            fileoutputstream = new FileOutputStream(fileName);
            dataoutputstream = new DataOutputStream(fileoutputstream);
            bytes = StsMath.floatsToUnsignedBytes254(attribute,aMin, aMax + onePrct);
            dataoutputstream.write(bytes);
            progressPanel.appendLine("Completed Output of crossline file: " + fileName);
            progressPanel.setValue(2);

            // Output Slices  -- Need to sort
            fileName = currentFile.substring(0, slashLoc+1) + "seis3d.bytes-trace." + currentFile.substring(slashLoc+1,dotLoc);
            progressPanel.appendLine("Outputting trace file: " + fileName);
            fileoutputstream = new FileOutputStream(fileName);
            dataoutputstream = new DataOutputStream(fileoutputstream);

            for(int n=0; n<numSlices; n++)
            {
                for(int i=0; i<numInlines*numXLines; i++)
                {
                    int offset = i*numSlices + n;
                    float attributeVal = attribute[offset];
                    byte byteVal = StsMath.floatToUnsignedByte254(attributeVal, aMin, aMax+onePrct);
                    dataoutputstream.write(byteVal);
                }
            }
            progressPanel.appendLine("Completed Output of trace file: " + fileName);
            progressPanel.setValue(3);
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
                aiFile.close();
                done = true;
            }
         }
      }
}
