package com.Sts.Actions.Wizards.HandVelocity;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.util.*;

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
public class StsHandVelocityWizard extends StsWizard
{

    StsCultureObjectSet2D cultureSet2D = null;
    float rowNumMin;
    float rowNumMax;
    float colNumMin;
    float colNumMax;

    protected StsHandVelocitySelect handVelSelect;
    protected StsHandVelocitySurveyDefinition handVelSurveyDef;

    private ArrayList<StsFile> handVelFilesList = new ArrayList<StsFile>();
    private StsSeismicBoundingBox[] handVelVolumes  = new StsSeismicBoundingBox[0];
    private StsSeismicBoundingBox[] selectedVolumes = new StsSeismicBoundingBox[0];

    public StsHandVelocityWizard(StsActionManager actionManager)
    {
        super(actionManager, 600, 600);
        addSteps
        (
            new StsWizardStep[]
            {
                handVelSelect    = new StsHandVelocitySelect(this),
                handVelSurveyDef = new StsHandVelocitySurveyDefinition(this)
            }
        );
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Load Hand Picked Velocity Profiles (HANDVEL)");
        dialog.getContentPane().setSize(500, 600);
        return super.start();
    }

    public boolean end()
    {
        if (success) model.setActionStatus(getClass().getName(), StsModel.STARTED);
        return super.end();
    }

    public void previous()
    {
//        if( currentStep == handVelSurveyDef)
//            removeFiles();
        
        gotoPreviousStep();
    }

    public void next()
    {
        if( currentStep == handVelSelect)
        {
            createSurveyVolumes();
            enableFinish();
        }
        gotoNextStep();
    }

    public void finish()
    {
        super.finish();
    }

    private boolean analyzeHandVelocity( StsAbstractFile file)
    {
        StringTokenizer st;
        String strToken;

        StsAsciiFile asciiFile = new StsAsciiFile(file);
        if(!asciiFile.openReadWithErrorMessage()) return false;

        try
        {
            String aline = asciiFile.readLine();
            while(aline != null)
            {
                if(aline.indexOf("HANDVEL") == 0)
                {
                    st = new StringTokenizer(aline, " ");
                    int nTokens = st.countTokens();
                    if(nTokens < 2)
                    {
                        new StsMessage(model.win3d,  StsMessage.WARNING,  "HandVel line: " + aline + " not sufficient to load.");
                        return false;
                    }
                    strToken = st.nextToken();
                    strToken = st.nextToken();
                    int inlineAndCrossline = Integer.valueOf(strToken);
                    float crossline = -1, inline = -1;
                    if(nTokens == 6) //extended Handvel format from Focus SWC 6/17/09
        			{
        				inline = Integer.parseInt(aline.substring(16,24).trim());
        				crossline = Integer.parseInt(aline.substring(25,32).trim());
        			}
                    if(strToken.length() > 4)
                    {
                        crossline = inlineAndCrossline%10000;
                        inline = (inlineAndCrossline - crossline)/10000;
                    }
                    if (nTokens == 2) {
                     // 2D Case!!!
                        inline = 1;
                        crossline = inlineAndCrossline;
                    }
                    else
                    { 
                        // 3D inline and xline separate case
                        inline = inlineAndCrossline;
                        strToken = st.nextToken();
                        crossline = Integer.valueOf(strToken);
                    }
                    rowNumMin = Math.min( inline, rowNumMin);
                    rowNumMax = Math.max( inline, rowNumMax);
                    colNumMin = Math.min( crossline, colNumMin);
                    colNumMax = Math.max( crossline, colNumMax);
                }
                aline = asciiFile.readLine();
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsHandVelocityWizard.analyzeHandVelocity() failed.", e, StsException.WARNING);
            return false;
        }
    }

    private void initializeDataRanges()
    {
        rowNumMin =  Float.MAX_VALUE;
	    rowNumMax = -Float.MAX_VALUE;
        colNumMin =  Float.MAX_VALUE;
	    colNumMax = -Float.MAX_VALUE;
    }

    private void createSurveyVolumes()
    {
        removeAllVolumes();

        int nSelected = handVelFilesList.size();
        for(int n = 0; n < nSelected; n++)
        {
            initializeDataRanges();
            StsFile handVelFile = (StsFile)handVelFilesList.get(n);
            if( !analyzeHandVelocity( handVelFile)) continue;

            StsSeismicBoundingBox surveyVolume =  new StsSeismicBoundingBox();
            surveyVolume.setRowNumMin( rowNumMin);
            surveyVolume.setRowNumMax( rowNumMax);
            surveyVolume.setColNumMin( colNumMin);
            surveyVolume.setColNumMax( colNumMax);
            surveyVolume.setNRows( (int)(rowNumMax - rowNumMin) + 1);
            surveyVolume.setNCols( (int)(colNumMax - colNumMin) + 1);
            String filename = handVelSelect.panel.filenameFilter.getFilenameStem( handVelFile.getFilename());
            surveyVolume.setName( filename);
            addSurveyVolume( surveyVolume);
        }
    }

    public void addFile( StsFile file) { handVelFilesList.add( file);}

    public void addSurveyVolume(StsSeismicBoundingBox surveyVolume)
    {
        handVelVolumes = (StsSeismicBoundingBox[]) StsMath.arrayAddElement(handVelVolumes, surveyVolume);
    }

    public void removeAllVolumes()
    {
        if(handVelVolumes == null) return;
        for(int i = 0; i<handVelVolumes.length; i++)
            handVelVolumes = (StsSeismicBoundingBox[]) StsMath.arrayDeleteElement(handVelVolumes, handVelVolumes[i]);

        handVelVolumes = new StsSeismicBoundingBox[0];
    }

    public void removeFile( StsFile file) { handVelFilesList.remove( file);}
    public void removeFiles() { handVelFilesList.clear();}


    public ArrayList<StsFile> getHandVelFilesList() { return handVelFilesList;}
    public StsSeismicBoundingBox[] getSelectedVolumes() { return selectedVolumes;}
    public StsSeismicBoundingBox[] getHandVelVolumes() { return handVelVolumes;}

    public StsSeismicBoundingBox getFirstSelectedVolume()
    {
        if(selectedVolumes == null || selectedVolumes.length == 0)  return null;
        return selectedVolumes[0];
    }

    public void setSelectedVolumes(StsSeismicBoundingBox[] selectedVolumes)
    {
        StsSeismicBoundingBox[] newSelectedVolumes = new StsSeismicBoundingBox[selectedVolumes.length];
        for(int i = 0; i < selectedVolumes.length; i++)
        {
            newSelectedVolumes[i] = selectedVolumes[i];
        }
        this.selectedVolumes = newSelectedVolumes;
    }

    public void setSelectedVolumes( int[] selectedIndices)
    {
        int selectedSize = Math.min( selectedIndices.length, handVelVolumes.length);
        StsSeismicBoundingBox[] newSelectedVolumes = new StsSeismicBoundingBox[selectedSize];
        for(int i = 0; i < selectedSize; i++)
            newSelectedVolumes[i] = handVelVolumes[selectedIndices[i]];
        this.selectedVolumes = newSelectedVolumes;
    }

    public ArrayList getHandVelVolumesList()
    {
        ArrayList handVelVolumesList = new ArrayList();
        for( int i = 0; i < handVelFilesList.size(); i++)
             handVelVolumesList.add( handVelVolumes[i]);
        return handVelVolumesList;
    }

//****wrw test  TODO this needs to be implemented for handVels
    public void moveFileToAvailableList() {}
    {
//        volumeSelect.addAvailableVolume( segyVolumes[i]);
    }
}