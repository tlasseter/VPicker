package com.Sts.Actions.Wizards.AncillaryData;

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

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.io.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsAncillaryDataSelectPanel extends StsJPanel implements StsFileTransferObjectFaceNew
{
    private StsFileTransferPanel selectionPanel;
    private AncillaryDataFilenameFilter filenameFilter = null;

    String[] suffix = null;
    private DefaultListModel availableVolsListModel = new DefaultListModel();

    private DefaultListModel selectedObjectsListModel = new DefaultListModel();
    private JTextField currentDirectoryField = new JTextField();

    private JList availableObjectsList = new JList();

    private StsAncillaryDataWizard wizard;
    private StsAncillaryDataSelect wizardStep;

    private String currentDirectory = null;
    private StsAbstractFileSet availableFileSet;
    private StsAbstractFile[] availableFiles;

    StsModel model = null;
    private String outputDirectory = null;
    private StsFile currentSelectedFile = null;

    private StsButton viewButton = null;
    private boolean archiveIt = true;

    public StsAncillaryDataSelectPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsAncillaryDataWizard)wizard;
        this.wizardStep = (StsAncillaryDataSelect)wizardStep;
        model = wizard.getModel();
        try
        {
            if (model != null)
            {
                StsProject project = model.getProject();
                currentDirectory = project.getRootDirString();
                outputDirectory = project.getDataFullDirString();
            }
            else
            {
                currentDirectory = System.getProperty("user.dirNo"); // standalone testing
                outputDirectory = currentDirectory;
            }

			suffix = StsAncillaryData.extensions[StsAncillaryDataClass.DOCUMENT];
			suffix = (String[])StsMath.arrayAddArray(suffix,StsAncillaryData.extensions[StsAncillaryDataClass.IMAGE]);
			suffix = (String[])StsMath.arrayAddArray(suffix,StsAncillaryData.extensions[StsAncillaryDataClass.MULTIMEDIA]);
			suffix = (String[])StsMath.arrayAddArray(suffix,StsAncillaryData.extensions[StsAncillaryDataClass.OTHER]);
			StsFilenameEndingFilter filenameFilter = new StsFilenameEndingFilter(suffix);
			selectionPanel = new StsFileTransferPanel(currentDirectory, filenameFilter, this, 350, 150);
            constructBeans();
            jbInit();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public boolean hasReloadButton()
    {
        return false;
    }
    public boolean hasDirectorySelection()
    {
        return true;
    }
    public void setReload(boolean reload)
    {}

    public boolean getReload()
    {
        return true;
	}
    public boolean hasArchiveItButton()
    {
        return true;
    }
    public boolean hasOverrideButton() { return false; }
    public void setOverrideFilter(boolean override) {}
    public boolean getOverrideFilter() { return false; }

    public void setArchiveIt(boolean archive)
    {
        archiveIt = archive;
    }

    public boolean getArchiveIt()
    {
        return archiveIt;
	}

    public void fileSelected(StsAbstractFile selectedFile)
    {
        if(selectedFile == null)
        {
            currentSelectedFile = null;
        }
        if (selectedFile != currentSelectedFile)
        {
            currentSelectedFile = (StsFile)selectedFile;
        }
    }
    
    public void availableFileSelected(StsAbstractFile selectedFile)
    {
    }
    
    public void addFiles(StsAbstractFile[] files)
    {
        for (int n = 0; n < files.length; n++)
            wizard.addFile((StsFile)files[n], outputDirectory);
        fileSelected(files[0]);
    }

    public void addAllFiles(StsFile[] files)
    {
        for (int n = 0; n < files.length; n++)
            wizard.addFile(files[n], outputDirectory);
        fileSelected(files[files.length-1]);
    }

    public void removeFiles(StsAbstractFile[] files)
    {
        for (int n = 0; n < files.length; n++)
            wizard.removeFile((StsFile)files[n]);
    }

    public void removeAllFiles()
    {
        wizard.removeFiles();
    }

    public void initialize()
    {
        filenameFilter = new AncillaryDataFilenameFilter(suffix);
        setAvailableFiles();
        availableObjectsList.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
    }

    private void setAvailableFiles()
    {
        availableFileSet = StsFileSet.constructor(currentDirectory, filenameFilter);
        availableFiles = availableFileSet.getFiles();
        int nFiles = availableFiles.length;
        String[] fileEndings = getFilenameEndings(availableFiles);
        int cnt = 0;
        availableVolsListModel.removeAllElements();
        if(fileEndings != null) {
          for(int i=0; i<fileEndings.length; i++) {
              if((StsSeismicVolume)model.getObjectWithName(StsSeismicVolume.class, fileEndings[i]) == null) {
                  availableVolsListModel.addElement(fileEndings[i]);
                  availableFiles[cnt] = availableFiles[i];
                  cnt++;
              }
          }
        }
    }

    private void constructBeans()
    {
        viewButton = new StsButton("View Selected File", "Review selected file contents", this, "showFile");
	}

    public void showFile()
    {
        String line;
        if((currentSelectedFile.getFilename().indexOf("txt") <= 0) && (currentSelectedFile.getFilename().indexOf("txt") <= 0))
        {
            boolean yesNo = StsYesNoDialog.questionValue(model.win3d, "Can only view text files, do you want to try anyway?");
            if(!yesNo)
                return;
        }
        StsAsciiFile selectedFile = new StsAsciiFile(currentSelectedFile);
        if(selectedFile == null)
        {
            new StsMessage(wizard.getModel().win3d, StsMessage.WARNING, "No file selected in list, must be in selected (left) list");
            return;
        }
        if(!selectedFile.openReadWithErrorMessage()) return;

        StsFileViewDialog dialog = new StsFileViewDialog(wizard.frame,"Hand Velocity File View", false);
        dialog.setVisible(true);
        dialog.setViewTitle("File - " + currentSelectedFile.getFilename());
        try {
            while (true)
            {
                line = selectedFile.readLine();
                if (line == null)
                    break;
                dialog.appendLine(line);
            }
        }
        catch (Exception ex)
        {
            StsMessageFiles.infoMessage("Unable to view selected file.");
        }
        finally
        {
            selectedFile.close();
        }
    }

    void jbInit() throws Exception
    {
        this.setLayout(new GridBagLayout());

        Color panelColor = Color.black;
        Font  panelFont = new java.awt.Font("Dialog", 1, 11);

        gbc.gridwidth = 2;
        gbc.fill = gbc.HORIZONTAL;
        addEndRow(selectionPanel);

        gbc.fill = gbc.NONE;
        addEndRow(viewButton);
    }

    public String getCurrentDirectory() { return currentDirectory; }

    public String[] getFilenameEndings(StsAbstractFile[] files)
    {
        int nFiles = files.length;
        String[] fileEndings = new String[nFiles];
        for(int n = 0; n < nFiles; n++)
            fileEndings[n] = filenameFilter.getFilenameStem(files[n].getFilename());
        return fileEndings;
    }

    final class AncillaryDataFilenameFilter implements FilenameFilter
     {
         String[] filter;
         int length;

         public AncillaryDataFilenameFilter(String[] filter)
         {
             this.filter = new String[filter.length];
             for (int i = 0; i < filter.length; i++)
                 this.filter[i] = filter[i].toLowerCase();
         }

         public boolean accept(File dir, String name)
         {
             for (int i = 0; i < filter.length; i++)
             {
                 if (name.toLowerCase().endsWith(filter[i]))
                     return true;
             }
             return false;
        }

        public String getFilenameStem(String filename)
        {
            int filterStart = 1;
            for (int i = 0; i < filter.length; i++)
            {
                if (filename.toLowerCase().indexOf(filter[i]) > 0)
                {
                    filterStart = filename.toLowerCase().indexOf(filter[i]);
                    break;
                }
            }
            return filename.substring(0, filterStart - 1);
		}
    }
}
