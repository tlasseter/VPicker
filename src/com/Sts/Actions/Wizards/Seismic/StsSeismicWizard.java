package com.Sts.Actions.Wizards.Seismic;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.util.*;

abstract public class StsSeismicWizard extends StsWizard implements Runnable
{
    protected StsPostStackSelect volumeSelect;
    protected StsPostStackFileFormat fileFormat;

    private boolean skipReanalyze = false;

    private double defaultScanPercent = 0.1;
    private transient boolean overrideHeader = false;
    private boolean overrideGeometry = false;
    protected transient StsSEGYFormat segyFormat;
    protected ArrayList segyVolumesList = new ArrayList();
    protected StsSeismicBoundingBox[] segyVolumesToProcess = new StsSeismicBoundingBox[0];
    protected StsSeismicBoundingBox[] selectedVolumes = new StsSeismicBoundingBox[0];
    protected String volumeName = null;

    public StsSeismicWizard(StsActionManager actionManager, int width, int height, byte segyType)
    {
        super(actionManager, width, height);
 //       setSegyFormat(StsSEGYFormat.constructor(model, segyType));
        addSteps();
        dialog.setDefaultCloseOperation(dialog.DO_NOTHING_ON_CLOSE);
    }


    abstract public StsPostStackAnalyzer getAnalyzer(StsProgressPanel panel, StsTablePanelNew volumeStatusPanel);
    abstract protected void addSteps();
    abstract public boolean start();
    abstract public byte getSegyFormatDataType();

    protected boolean start(String title)
    {
        System.runFinalization();
        System.gc();
        disableFinish();
        dialog.setTitle(title);
        return super.start();
    }

    public void setWizardSize(int w, int h)
    {
        dialog.getContentPane().setSize(w, h);
        dialog.validate();
    }

    abstract public void previous();

    abstract public void next();

    protected void saveSEGYFormats()
    {
        if(!segyFormat.isFormatChanged()) return;
        StsSEGYFormatSaveDialog dialog = new StsSEGYFormatSaveDialog(model.win3d, new StsSEGYFormat[] { segyFormat } );
        dialog.setVisible(true);
    }

    public boolean getIsLittleEndian()
     {
         return segyFormat.getIsLittleEndian();
     }

    public StsSeismicBoundingBox[] getSegyVolumes()
    {
        if( segyVolumesList.size() <= 0) return new StsSeismicBoundingBox[0];

        StsSeismicBoundingBox[] volumes = new StsSeismicBoundingBox[0];
        for( int i = 0; i < segyVolumesList.size(); i++)
            volumes = (StsSeismicBoundingBox[]) StsMath.arrayAddElement(volumes, segyVolumesList.get(i));
        return volumes;
    }

    public StsSeismicBoundingBox[] getVolumes()
    {
        return getSegyVolumes();
    }

    public ArrayList getSegyVolumesList()
    {
        return segyVolumesList;
    }

    public StsSeismicBoundingBox[] getSegyVolumesToProcess()
    {
        return segyVolumesToProcess;
    }

    public ArrayList getSegyVolumesToProcessList()
    {
        ArrayList segyVolumesToProcessList = new ArrayList();
        for( int i = 0; i < segyVolumesToProcess.length; i++)
             segyVolumesToProcessList.add( segyVolumesToProcess[i]);
        return segyVolumesToProcessList;
    }

    abstract public boolean isStatusOk( StsSeismicBoundingBox segyVolume);

    public void addSegyVolumesToProcessList()
    {
        StsSeismicBoundingBox[] segyVolumes = getSegyVolumes();
        for( int i = 0; i < segyVolumes.length; i++)
        {
            if( isStatusOk(segyVolumes[i]))
                segyVolumesToProcess = (StsSeismicBoundingBox[]) StsMath.arrayAddElement( segyVolumesToProcess, segyVolumes[i]);
            else
                volumeSelect.addAvailableVolume( segyVolumes[i]);
        }

        selectedVolumes = null;
        segyVolumesList.clear();
        volumeSelect.clearFiles();
    }

    public boolean end()
    {
        model.setActionStatus(getClass().getName(), StsModel.STARTED);
        return super.end();
    }

    public void removeSegyVolumeToProcess(StsSeismicBoundingBox volume)
    {
        segyVolumesToProcess = (StsSeismicBoundingBox[])StsMath.arrayDeleteElement(segyVolumesToProcess, volume);
        volumeSelect.addAvailableVolume(volume);
    }

    public StsSeismicBoundingBox getSegyVolume( int index)
    {
        StsSeismicBoundingBox volume = (StsSeismicBoundingBox)segyVolumesList.get(index);
        return volume;
    }

    public StsSeismicBoundingBox getSegyVolume(String filename)
    {
        StsSeismicBoundingBox[] volumes = getSegyVolumes();
        for(int i = 0; i < volumes.length; i++)
        {
            if(volumes[i].getSegyFilename().equals(filename))
                return volumes[i];
        }
        return null;
    }

    public void setZDomainString(String domain)
    {
        if(segyFormat == null) return;
        
        segyFormat.setZDomainString(domain);
        if(!segyVolumesSEGYFormatAvailable()) return;
        StsSeismicBoundingBox[] volumes = getSelectedVolumes();
        for(int i = 0; i < volumes.length; i++)
            volumes[i].setZDomainString(domain);
    }

    private boolean segyVolumesSEGYFormatAvailable()
    {
        if(segyVolumesList.size() == 0) return false;

        StsSeismicBoundingBox segyVolume = (StsSeismicBoundingBox)segyVolumesList.get(0);
        if( segyVolume.getSegyFormat() == null) return false;

        return true;
    }

    public String getZDomainString()
    {
        if(segyFormat == null) return model.getProject().getZDomainString();
        return segyFormat.getZDomainString();
    }

    public void setHorzUnits(String units)
    {
        if(segyFormat == null) return;

        segyFormat.setHUnitString(units);
        float hScalar = model.getProject().getXyScalar(units);
        StsSeismicBoundingBox[] volumes = getSelectedVolumes();
        for(int i = 0; i < volumes.length; i++)
            volumes[i].setHorizontalScalar(hScalar);
    }

    public String getHorzUnits()
    {
        if(segyFormat == null) return model.getProject().getXyUnitString();
        return segyFormat.hUnitString;
    }

    public void setVertUnits(String units)
    {
        if(segyFormat == null) return;

        float vScalar;
        if(segyFormat.zDomain == StsParameters.TD_TIME)
        {
            segyFormat.setTUnitString(units);
            vScalar = model.getProject().getTimeScalar(units);
        }
        else
        {
            segyFormat.setZUnitString(units);
            vScalar = model.getProject().getDepthScalar(units);
        }

        StsSeismicBoundingBox[] volumes = getSelectedVolumes();
        for(int i = 0; i < volumes.length; i++)
            volumes[i].setVerticalScalar(vScalar);
    }

    public String getVertUnits()
    {
        if(segyFormat == null) return model.getProject().getDepthUnitString();
        if(segyFormat.zDomain == StsParameters.TD_TIME)
            return segyFormat.tUnitString;
        else
            return segyFormat.zUnitString;
    }

    public void saveSelectedAttributes(StsSeismicBoundingBox volume) {}

    public StsSegyInformationPanel getInfoPanel()
    {
        return null;
        //   return infoPanel;
    }

    public void clearSegyVolumesList()
    {
        segyVolumesList.clear();
    }

    public StsSeismicBoundingBox addSegyVolume(StsFile file, String outputDirectory)
    {
        StsSeismicBoundingBox volume = createNewStsSegyVolume(this, file, outputDirectory, segyFormat);
        if(volume != null) segyVolumesList.add(volume);
        return volume;
    }

    public void checkForOutputFiles()
    {
        int nFiles = segyVolumesList.size();
        StsSeismicBoundingBox[] volumesWithOutput = new StsSeismicBoundingBox[nFiles];
        int nFilesWithOutput = 0;
        for(int n = 0; n < nFiles; n++)
        {
            StsSeismicBoundingBox volume = (StsSeismicBoundingBox)segyVolumesList.get(n);
            if(volume.hasOutputFiles()) volumesWithOutput[nFilesWithOutput++] = volume;
        }
        if(nFilesWithOutput == 0) return;
        StsProgressTextPanel textPanel = StsProgressTextPanel.constructor(nFilesWithOutput+3, 30);
        textPanel.appendLine("The following volumes have output:");
        for(int n = 0; n < nFilesWithOutput; n++)
            textPanel.appendLine("    " + volumesWithOutput[n].getName());
        textPanel.appendLine("Delete any of the files you don't want to process.");
        new StsOkDialog(frame, textPanel, "Already processed files", true);
    }

    protected StsSeismicBoundingBox createNewStsSegyVolume(StsSeismicWizard seismicWizard, StsFile file, String outputDirectory, StsSEGYFormat format)
    {
        return StsSegyVolume.constructor(seismicWizard, file, outputDirectory, frame, format);
    }

    public void removeSegyVolume(StsSeismicBoundingBox volume)
    {
        segyVolumesList.remove( volume);
    }

    public void moveSegyVolumeToAvailableList( StsSeismicBoundingBox volume)
    {
        volumeSelect.moveToAvailableList(volume);
    }

    public StsSEGYFormat getSegyFormat()
    {
        return segyFormat;
    }

    public void resetSelectedTraceStatus()
    {
      if( selectedVolumes == null) return;
      for( int i = 0; i < selectedVolumes.length; i++)
          selectedVolumes[i].resetTraceStatus();
    }

    public void analyzeHeaders()
    {
        if(!skipReanalyze) fileFormat.reanalyzeHeaders();
    }

    public void analyzeTraces()
    {
        if(!skipReanalyze) fileFormat.reanalyzeTraces();
    }

    public StsSegyData[] getSelectedSegyDatasets()
    {
        return StsSegyData.getSegyDatasets(getSelectedVolumes());
    }

    public StsSegyData[] getVolumeSegyDatasets()
    {
        return StsSegyData.getSegyDatasets(getVolumes());
    }

    public StsSeismicBoundingBox getFirstSegyVolume()
    {
        if(segyVolumesList.size() == 0)
            return null;
        else
            return (StsSeismicBoundingBox)segyVolumesList.get(0);
    }

    public boolean resetOverrideValues()
    {
        boolean reset = false;
        if(!segyVolumesSEGYFormatAvailable())
            return reset;

        StsSeismicBoundingBox segyVolume = getFirstSegyVolume();
        if( segyVolume == null) return reset;

        int oldOverrideNSamples = segyFormat.getOverrideNSamples();
        int newOverrideNSamples = segyFormat.getNSamples(this.getIsLittleEndian(), segyVolume.getBinaryHeader());
        if(oldOverrideNSamples != newOverrideNSamples)
        {
            reset = true;
            segyFormat.setOverrideNSamples(newOverrideNSamples);
        }

        float oldOverrideSampleSpacing = segyFormat.getOverrideSampleSpacing();
        float newOverrideSampleSpacing = segyFormat.getSampleSpacing(this.getIsLittleEndian(), segyVolume.getBinaryHeader());
        if(oldOverrideSampleSpacing != newOverrideSampleSpacing)
        {
            reset = true;
            segyFormat.setOverrideSampleSpacing(newOverrideSampleSpacing);
        }

        if(reset)
           analyzeHeaders();

        return reset;
    }

    public boolean getOverrideHeader()
    {
        return overrideHeader;
    }

    public void setOverrideHeader(boolean override)
    {
        overrideHeader = override;
    }

    public void setSegyFormat(StsSEGYFormat segyFormat)
    {
        this.segyFormat = segyFormat;
    }

    public void setVolumesSegyFormat()
    {
        StsSeismicBoundingBox[] volumes = getSegyVolumes();
        for(int i = 0; i < volumes.length; i++)
            volumes[i].setSegyFormat(segyFormat);
    }

    public void setIsLittleEndian(boolean isLittleEndian)
    {
        segyFormat.setIsLittleEndian(isLittleEndian);
        if(!segyVolumesSEGYFormatAvailable()) return;
        StsSeismicBoundingBox[] volumes = getSegyVolumes();
        for(int i = 0; i < volumes.length; i++)
            volumes[i].setIsLittleEndian(isLittleEndian);
        analyzeHeaders();
    }

    public float calcDiskRequired()
    {
        return StsSeismicBoundingBox.calcDiskRequired(getSegyVolumes());
    }

    public void setVolumeName(String name)
    {
        volumeName = name;
    }

    public String getVolumeName()
    {
        return volumeName;
    }

    public double getScanPercent()
    {
        return defaultScanPercent;
    }

    public void setScanPercent(double scanPercent)
    {
        defaultScanPercent = scanPercent;
    }

    public double getDefaultScanPercent()
    {
        return defaultScanPercent;
    }

    public void setDefaultScanPercent(double defaultScanPercent)
    {
        this.defaultScanPercent = defaultScanPercent;
    }

    public void setSkipReanalyzeTraces(boolean value)
    {
        skipReanalyze = value;
    }

    public void setSelectedVolumes(StsSeismicBoundingBox[] selectedVolumes)
    {
        StsSeismicBoundingBox[] newSelectedVolumes = new StsSeismicBoundingBox[selectedVolumes.length];
        for(int i = 0; i < selectedVolumes.length; i++)
        {
            newSelectedVolumes[i] = selectedVolumes[i];
        }
        this.selectedVolumes = newSelectedVolumes;
//        updateFileFormatPanel();
    }

    public void setSelectedVolumes(int[] selectedIndices)
    {
        int selectedSize = Math.min( selectedIndices.length, segyVolumesList.size());
        StsSeismicBoundingBox[] newSelectedVolumes = new StsSeismicBoundingBox[selectedSize];
        for(int i = 0; i < selectedSize; i++)
            newSelectedVolumes[i] = (StsSeismicBoundingBox)segyVolumesList.get(selectedIndices[i]);
        this.selectedVolumes = newSelectedVolumes;
//        updateFileFormatPanel();
    }

    public void setSelectedVolumesToAll()
    {
        int selectedSize = segyVolumesList.size();
        StsSeismicBoundingBox[] newSelectedVolumes = new StsSeismicBoundingBox[selectedSize];
        for(int i = 0; i < selectedSize; i++)
            newSelectedVolumes[i] = (StsSeismicBoundingBox)segyVolumesList.get(i);
        this.selectedVolumes = newSelectedVolumes;
//        updateFileFormatPanel();
    }

    public StsSeismicBoundingBox[] getSelectedVolumes()
    {
        if(selectedVolumes == null) return new StsSeismicBoundingBox[0];
        return selectedVolumes;
    }

    public StsSeismicBoundingBox getFirstSelectedVolume()
    {
        if(selectedVolumes.length == 0) return null;
        return selectedVolumes[0];
    }

    public StsSeismicBoundingBox getSurveyDefinitionBoundingBox()
    {
        StsException.systemError(this, "getSurveyDefinitionBoundingBox", "Should be overridden in subclass.");
        return null;
    }

    public void updateFileFormatPanel()
    {
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                fileFormat.updateFileFormatPanel();
            }
        };
        StsToolkit.runWaitOnEventThread(runnable);
    }

    public void setOverrideGeometry(boolean override)
    {
        if(override == overrideGeometry) return;
        if(!StsRotatedGridBoundingBox.allAreCongruent( getSegyVolumes()))
        {
            new StsMessage(dialog, StsMessage.WARNING, "Selected volumes are not congruent,\n so override geometry cannot be applied.\n" +
                "Rerun wizard with a single volume or a group of congruent volumes.");
            return;
        }
        overrideGeometry = override;
        if(overrideGeometry) enableNext();
 //       disableNext();
 //       analyzeGrid();
    }

    public boolean getOverrideGeometry() { return overrideGeometry; }

    public boolean cropEnabled() { return true;}
}
