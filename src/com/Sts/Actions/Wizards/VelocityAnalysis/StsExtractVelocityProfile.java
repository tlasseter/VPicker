package com.Sts.Actions.Wizards.VelocityAnalysis;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.util.*;

public class StsExtractVelocityProfile extends StsWizardStep implements Runnable
{
    public StsProgressPanel panel;
    private StsHeaderPanel header;
    private StsVelocityAnalysisWizard wizard;
    private StsPreStackVelocityModel velocityModel = null;
    private StsPreStackVelocityModel inputVelocityModel = null;
	private boolean canceled = false;
    private float velocityScaleMultiplier = 1.0f;

    public StsExtractVelocityProfile(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsVelocityAnalysisWizard)wizard;
        constructPanel();
    }

    public void constructPanel()
    {
        panel = StsProgressPanel.constructorWithCancelButton(5, 50);
        header = new StsHeaderPanel();
        setPanels( panel, header);
        header.setTitle("Extract Velocity Profile");
        header.setSubtitle("Process Velocity Profile(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VelocityAnalysis");
        header.setInfoText(wizardDialog,"(1) Once loading is complete, press the Finish Button to dismiss the screen");
    }

    public boolean start()
    {
        wizard.dialog.setTitle("Process Velocity Profile(s)");
        run();
        return true;
    }

    private boolean processVolumeProfile( StsPreStackVelocityModel velocityModel)
    {
        if (velocityModel == null) return false;

        StsSeismicVolume velocityVolume = (StsSeismicVolume)wizard.getInputVelocityModelOrVolume();
      
        if( velocityVolume == null) return false;

        StsPreStackLineSet currentLineSet = wizard.getCurrentLineSet();
        if (currentLineSet == null) return false;

        int numVelPts = wizard.getNumVelPts();
        int nSamples = velocityVolume.getNSlices();
        int nVelPts = Math.min( nSamples, numVelPts);
        if( nVelPts < 1) return false;
        int sampleInc = Math.max((nSamples - 1)/nVelPts, 1);
        float zInc = velocityVolume.getZInc();
        velocityScaleMultiplier = (float)wizard.getVelScaleMultiplier();

        StsObjectList velocityProfileList = new StsObjectList(0);

        Iterator rowColIterator = currentLineSet.getRowColIterator(1);
        int nRows = velocityModel.getNRows();
        panel.setMaximum(nRows);

        panel.appendLine("Extracting profiles for model " + velocityModel.getName());

        int currentRow = -1;
        while(rowColIterator.hasNext())
        {
            int[] rowCol = (int[])rowColIterator.next();
            if( rowCol == null) return false;
            int lineSetRow = rowCol[0];
            int lineSetCol = rowCol[1];
            float lineSetRowNum = currentLineSet.getRowNumFromRow(lineSetRow);
            float lineSetColNum = currentLineSet.getColNumFromCol(lineSetCol);
            int volumeRow = velocityVolume.getRowFromRowNum(lineSetRowNum);
            int volumeCol = velocityVolume.getColFromColNum(lineSetColNum);
            volumeRow = Math.min( volumeRow, velocityVolume.getNRows() - 1);
            volumeCol = Math.min( volumeCol, velocityVolume.getNCols() - 1);
            if(!velocityVolume.isInsideRowColRange(volumeRow, volumeCol)) continue;
            float[] traceVelocity = velocityVolume.getTraceValues(volumeRow, volumeCol);
             if(traceVelocity == null) continue;

            StsVelocityProfile velocityProfile = new StsVelocityProfile(lineSetRow, lineSetCol);
            velocityProfile.initializeMutes(currentLineSet);
            for( int i = 0; i <= nVelPts; i++)
            {
                StsPoint point = new StsPoint(velocityScaleMultiplier*traceVelocity[i*sampleInc], i*sampleInc*zInc);
                velocityProfile.addProfilePoint(point);
            }
            velocityProfileList.add( velocityProfile);
            if(lineSetRow != currentRow)
            {
                panel.appendLine("Extracted profiles through row "+ lineSetRow + " of " + nRows);
                panel.setDescription(" " + velocityModel.getName() + ":" + lineSetRow);
                panel.setValue(lineSetRow);
                currentRow = lineSetRow;
            }
        }
        velocityModel.setVelocityProfiles( velocityProfileList.convertListToRefList(model, "velocityProfiles", velocityModel));
        if(!velocityModel.checkInitializeInterpolation()) return false;
        velocityModel.initializeInterpolateVelocityProfiles();

        panel.setValue( nRows);
        model.viewObjectRepaint(this, velocityModel);
        return true;
    }

    private int getAssocColIndex( StsSeismicLineSet assocLineSet, int assocLineIndex, int cdp)
    {
        StsSeismicLine2d line = assocLineSet.lines[assocLineIndex];
        int nCols = line.getNCols();
        if( nCols <= 0) return -1;

        for( int i=0; i<nCols; i++)
        {
            if( line.cdp[i] == cdp)
                return i;
        }
        return -1;
    }

    private int getAssocLineIndex( StsSeismicLineSet assocLineSet, String lineName)
    {
        int nRows = assocLineSet.getNRows();
        if( nRows <= 0) return -1;

        for( int i=0; i<nRows; i++)
        {
            StsSeismicLine seismicLine = assocLineSet.lines[i];
            if( seismicLine.getAssocLineName().compareTo(lineName) == 0)
                return i;
        }
        return -1;
    }

    private boolean processLineSetProfile( StsPreStackVelocityModel velocityModel)
    {
        if (velocityModel == null) return false;

        StsSeismicLineSet assocLineSet = (StsSeismicLineSet)wizard.getInputVelocityModelOrVolume();
        if( assocLineSet == null) return false;

        StsPreStackLineSet currentLineSet = wizard.getCurrentLineSet();
        if (currentLineSet == null) return false;

        int nLines = currentLineSet.getNRows();
        if( nLines < 1) return false;

        int numVelPts = wizard.getNumVelPts();
        velocityScaleMultiplier = (float)wizard.getVelScaleMultiplier();

        int[] assocLineIndices = new int[nLines];
        for( int i=0; i<nLines; i++)
        {
            String lineName = currentLineSet.lines[i].getName();
            assocLineIndices[i] = getAssocLineIndex( assocLineSet, lineName);
        }

        StsObjectList velocityProfileList = new StsObjectList(0);

        Iterator rowColIterator = currentLineSet.getRowColIterator(1);
        int nRows = velocityModel.getNRows();
        panel.setMaximum(nRows);

        panel.appendLine("\n");
        panel.appendLine("Extracting profiles for model " + velocityModel.getName());

        int nSamples = 0;
        int nVelPts = 0;
        int sampleInc = 0;
        float assocLineZInc = 0.0f;
        float[] traceVelocity = null;
        int currentRow = -1;
        while(rowColIterator.hasNext())
        {
            int[] rowCol = (int[])rowColIterator.next();
            if( rowCol == null) return false;
            int lineSetRow = rowCol[0];
            int lineSetCol = rowCol[1];

            int assocLineIndex = assocLineIndices[lineSetRow];
            if( assocLineIndex < 0) continue;

            StsPreStackLine currentLine = currentLineSet.lines[lineSetRow];
            StsSeismicLine2d assocLine  = assocLineSet.lines[assocLineIndex];

            if( lineSetRow != currentRow)
            {
                traceVelocity = assocLine.readFloatData();
                if(traceVelocity == null) continue;

                nSamples = assocLine.getNSlices();
                nVelPts = Math.min( nSamples, numVelPts);
                if( nVelPts < 1) continue;
                assocLineZInc = assocLineSet.getZInc();
                sampleInc = Math.max((nSamples - 1)/nVelPts, 1);
            }

            int currentCdp = currentLine.cdp[lineSetCol];
            int assocCol = getAssocColIndex( assocLineSet, assocLineIndex, currentCdp);
            if( assocCol < 0) continue;

            StsVelocityProfile velocityProfile = new StsVelocityProfile(lineSetRow, lineSetCol);
            velocityProfile.initializeMutes(currentLineSet);
            for( int i=0; i<=nVelPts; i++)
            {
                int velocityOffset = assocCol*nSamples + i*sampleInc;
                StsPoint point = new StsPoint(velocityScaleMultiplier*traceVelocity[velocityOffset], i*sampleInc*assocLineZInc);
                velocityProfile.addProfilePoint(point);
            }
            velocityProfileList.add( velocityProfile);
            if(lineSetRow != currentRow)
            {
                panel.appendLine("Extracted profiles through row "+ lineSetRow + " of " + nRows);
                panel.setDescription(" " + velocityModel.getName() + ":" + lineSetRow);
                panel.setValue(lineSetRow);
                ((StsPreStackLine2d)currentLine).textureChanged();

                currentRow = lineSetRow;
            }
        }
        velocityModel.setVelocityProfiles( velocityProfileList.convertListToRefList(model, "velocityProfiles", velocityModel));
        if(!velocityModel.checkInitializeInterpolation()) return false;
        velocityModel.initializeInterpolateVelocityProfiles();
        
        panel.setValue( nRows);
        model.viewObjectRepaint(this, velocityModel);

        return true;
    }

    private StsObjectList loadHandVelocity( StsAsciiFile asciiFile, StsPreStackVelocityModel velocityModel,
                                            boolean is3d, StsPreStackLineSet3d velocityVolume, StsPreStackLine assocLine,
                                            StsObjectList velocityProfileList)
    {
        StringTokenizer st;
        String strToken;
        boolean startHandvel = false;
        StsVelocityProfile velocityProfile = null;

        StsPreStackLineSet currentLineSet = wizard.getCurrentLineSet();
        if (currentLineSet == null) return velocityProfileList;

        int row = -1;
        int col = -1;
        int lineSetRow = -1;
        int currentRow = -1;
        try
        {
            String aline = asciiFile.readLine();
            while( aline != null)
            {
                if( aline.indexOf("HANDVEL") == 0)
                {
                    // Add previous profile to model
                    if( velocityProfile != null)
                        velocityProfileList.add( velocityProfile);

                    // Extract the inline and crossline
                    st = new StringTokenizer( aline, " ");
                    int nTokens = st.countTokens();
                    if(nTokens < 2) continue;
                    strToken = st.nextToken();
                    strToken = st.nextToken();
                    int inlineAndCrossline = Integer.valueOf(strToken).intValue();

                    if( is3d)
                    {
                        float crossline = -1, inline = -1;
                        if(nTokens == 6)  // extended Handvel format from Focus - SWC 6/18/09
                        {
                        	inline = Integer.parseInt(aline.substring(17,24).trim());
                        	crossline = Integer.parseInt(aline.substring(25,32).trim());
                        }
                        else if(strToken.length() > 4)
                        {
                            crossline = inlineAndCrossline%10000;
                            inline = (inlineAndCrossline - col)/10000;
                        }
                        else
                        {
                            if(nTokens < 4) continue;
                            inline = inlineAndCrossline;
                            strToken = st.nextToken();
                            crossline = Integer.valueOf(strToken);
                        }
                        row = currentLineSet.getRowFromRowNum(inline);
                        col = currentLineSet.getColFromColNum(crossline);
                        panel.appendLine("Processing 3d handVel profile at inline: " + inline + " crossline: " + crossline);
                    }
                    else
                    {
                        float crossline = inlineAndCrossline%10000;
                        row = assocLine.lineIndex;
                        col = StsSegyLine2d.getAssocCdpIndex( assocLine.cdp,  (int)crossline);
                        panel.appendLine("Processing 2d handVel profile at cdp: " + crossline);
                    }
                    lineSetRow = row;

                    if(velocityModel.isInsideRowColRange(row, col))
                    {
                        velocityProfile = velocityModel.getVelocityProfile( row, col);
                        if( velocityProfile != null)
                            velocityProfile.clearProfilePoints();
                        else
                        {
                            velocityProfile = new StsVelocityProfile( row, col);
                            velocityProfile.initializeMutes(currentLineSet);
                            startHandvel = true;
                        }
                    }
                    else
                    {
                        startHandvel = false;
                        panel.appendLine("   handVel is outside volume range " + velocityModel.getRowColNumLabel());    
                    }

                    if( startHandvel)
                        aline = asciiFile.readLine();
                }
                if( startHandvel)
                {
                    st = new StringTokenizer( aline, " ");
                    while( st.hasMoreTokens())
                    {
                        strToken = st.nextToken();
                        float y = Float.valueOf(strToken).floatValue();
                        strToken = st.nextToken();
                        float velocity = Float.valueOf(strToken).floatValue()*velocityScaleMultiplier;
                        velocityProfile.addProfilePoint( velocity, y);
                    }

                    if( is3d && lineSetRow != currentRow)
                    {
                        panel.appendLine( "Extracted profiles through row "+ lineSetRow + " of " + velocityModel.getNRows());
                        panel.setDescription( " " + velocityModel.getName() + ":" + lineSetRow);
                        panel.setValue( lineSetRow);
                        currentRow = lineSetRow;
                    }
                }
                aline = asciiFile.readLine();
            }
            if( velocityProfile != null)
                velocityProfileList.add( velocityProfile);
            return velocityProfileList;
        }
        catch(Exception e)
        {
            String fileLine = asciiFile.getLine();
            StsException.systemError("Failed to read line for file line:\n" + fileLine);
            return velocityProfileList;
        }
    }

    private StsFile[] addHandVelFile( String directory, String handVelFilename, StsFile[] handVelFiles,
                                      StsFilenameEndingFilter filenameFilter)
    {
        String filename = handVelFilename + ".handvel";
        StsFile newFile = StsFile.constructor( directory, filename, filenameFilter);
        handVelFiles = (StsFile[] )StsMath.arrayAddElement( handVelFiles, newFile);
        return handVelFiles;
    }

    private StsFile[] getHandVelFiles( boolean is3d)
    {
        /*
        String[] suffix = new String[] {"handvel"};
        StsFilenameEndingFilter filenameFilter = new StsFilenameEndingFilter(suffix);
         */
        StsFile[] handVelFiles = new StsFile[0];
        if( is3d)
        {
            Object[] volumes = (Object[] )wizard.getModel().getTrimmedList( StsPreStackLineSet3d.class);
            for( int vol=0; vol<volumes.length; vol++)
            {
                StsPreStackLineSet3d seismicVolume = (StsPreStackLineSet3d)volumes[vol];
                if( seismicVolume.getHandVelName() != null)
                    handVelFiles = seismicVolume.getHandVelFiles();
            }
        }
        else
        {
            Object[] preStackLineSet = (Object[] )wizard.getModel().getTrimmedList( StsPreStackLineSet2d.class);
            for( int pset=0; pset<preStackLineSet.length; pset++)
            {
                StsPreStackLineSet pLineSet = (StsPreStackLineSet2d)preStackLineSet[pset];
                for( int pline=0; pline<pLineSet.lines.length; pline++)
                {
                    StsPreStackLine preStackLine = pLineSet.lines[pline];
                    if( preStackLine.getHandVelName() != null)
                        handVelFiles = preStackLine.getHandVelFiles();
                }
            }
        }
        return handVelFiles;
    }

    private StsPreStackLine[] getLineForHandVel( StsPreStackLineSet lineSet, StsFile[] handVelFiles)
    {
        StsPreStackLine[] lineForHandVel = new StsPreStackLine[handVelFiles.length];
        if( handVelFiles.length <= 0)
            return lineForHandVel;

        int nLines = lineSet.getNRows();
        if( nLines < 1) return lineForHandVel;

        for( int i=0; i<nLines; i++)
        {
            StsPreStackLine preStackLine = lineSet.lines[i];
            if( preStackLine.getHandVelName() == null) continue;
            for( int j=0; j<handVelFiles.length; j++)
            {
                if( handVelFiles[j].getFilenameStem().compareTo(preStackLine.getHandVelName()) == 0)
                {
                    lineForHandVel[j] = preStackLine;
                    break;
                }
            }
        }
        return lineForHandVel;
    }

    private boolean processHandVelProfile( StsPreStackVelocityModel velocityModel)
    {
        if (velocityModel == null) return false;

        boolean is3d = false;
        if( wizard.getCurrentLineSet() instanceof StsPreStackLineSet3d)
            is3d = true;

        StsFile[] handVelFiles;
        handVelFiles = getHandVelFiles( is3d);
        int nFiles = handVelFiles.length;
        if( nFiles <= 0) return true;

        StsPreStackLineSet3d velocityVolume = null;
        StsPreStackLine[] assocLines = new StsPreStackLine[nFiles];
        int maxCount;
        if( is3d)
        {
            velocityVolume = (StsPreStackLineSet3d)wizard.getCurrentLineSet();
            if( velocityVolume == null) return false;
            maxCount = velocityModel.getNRows();
        }
        else
        {
            StsPreStackLineSet currentLineSet = wizard.getCurrentLineSet();
            if (currentLineSet == null) return false;
            assocLines = getLineForHandVel( currentLineSet, handVelFiles);
            maxCount = nFiles;
        }

        if( !velocityModel.checkInitializeInterpolation()) return false;

        StsObjectList velocityProfileList = new StsObjectList(0);
        velocityScaleMultiplier = (float)wizard.getVelScaleMultiplier();

        panel.setMaximum( maxCount);
        panel.appendLine( "\n");
        panel.appendLine( "Extracting  handVel profiles ");

        for( int file=0; file<nFiles; file++)
        {
            StsAsciiFile asciiFile = new StsAsciiFile( handVelFiles[file] );
            if(!asciiFile.openReadWithErrorMessage())
            {
                panel.appendLine( "Error reading file " + handVelFiles[file].getFilename());
                continue;
            }

            velocityProfileList = loadHandVelocity( asciiFile, velocityModel, is3d, velocityVolume, assocLines[file],
                                                    velocityProfileList);
            if( !is3d)
            {
                panel.appendLine( "\n");
                panel.appendLine( "Extracted profiles for handVel file " + handVelFiles[file].getFilename());
                panel.setDescription("Completed handVel file " + handVelFiles[file].getFilename());
                panel.setValue( file + 1);
                ((StsPreStackLine2d)assocLines[file]).textureChanged();
            }
            asciiFile.close();
        }

        velocityModel.setVelocityProfiles( velocityProfileList.convertListToRefList(model, "velocityProfiles", velocityModel));
        velocityModel.initializeInterpolateVelocityProfiles();
        model.viewObjectRepaint(this, velocityModel);

        return true;
   }

    public boolean doProcessProfile( StsPreStackVelocityModel velocityModel)
    {
        Object modelType = wizard.getInputVelocityModelOrVolume();
        if( modelType instanceof StsSeismicLineSet)
            return processLineSetProfile( velocityModel);
        else if( modelType instanceof StsPreStackVelocityModel )
            return processPreStackProfile( velocityModel );
        else if( modelType instanceof StsSeismic )
            return processVolumeProfile( velocityModel);
        else if( modelType == StsPreStackLineSet.HAND_VEL)
            return processHandVelProfile( velocityModel);

        return false;
    }

    /**
     * Extracts velocity profiles from previous StsPreStackVelocityModel
     * 
     * right now it only works if output is within bounds of input
     * 
     * @param output
     * @return true for success, false for failure
     */
    private boolean processPreStackProfile(StsPreStackVelocityModel output)
    {
        if (output == null) return false;
        if (!output.checkInitializeInterpolation()) return false;
        if (velocityModel == output) return true; //keeps this from running twice
        
        StsPreStackVelocityModel input = (StsPreStackVelocityModel)wizard.getInputVelocityModelOrVolume();
      
        if( input == null) return false;

        StsPreStackLineSet currentLineSet = wizard.getCurrentLineSet();
        if (currentLineSet == null) return false;
        output.lineSet = currentLineSet;
        
        Object[] inputProfiles = input.getVelocityProfiles().getArrayList();
        for (int i=0; i< inputProfiles.length; i++)
        {
            StsVelocityProfile inputProfile = (StsVelocityProfile)inputProfiles[i];
            if (inputProfile != null && inputProfile.getProfilePoints().length > 1)
            {
                //... Convert column of input profile to what it will be in output velocity model
                int inputCol = inputProfile.col;
                int inputRow = inputProfile.row;
                float[] inputXY = input.getXYCoors(inputRow, inputCol);
                double[] inputAbsXY = input.getAbsoluteXY(inputXY[0], inputXY[1]);
                float[] outputXY = output.getRotatedRelativeXYFromUnrotatedAbsoluteXY(inputAbsXY[0], inputAbsXY[1]);
                int[] outputRowCol = output.lineSet.getRowColFromCoors(outputXY[0], outputXY[1]);
                if ( outputRowCol[0] >= 0 && outputRowCol[1] >= 0 )
                {
                    inputProfile.row = outputRowCol[0];
                    inputProfile.col = outputRowCol[1];
                    output.addProfile(inputProfile);
                }
                else
                {
                    StsMessage.printMessage("StsExtractVelocityProfile.processPreStackProfile : input profile outside bounds of output model: " + input.getName() + 
                            " row: " + inputRow + " col: " + inputCol);
                }
            }
        }
        if(!output.checkInitializeInterpolation()) return false;
        output.initializeInterpolateVelocityProfiles();
        
        currentLineSet.setVelocityModel(output);
        velocityModel = output;
        currentLineSet.superGathers[0].reinitialize();
        
        panel.setValue( currentLineSet.nRows );
        model.viewObjectChanged(this, output);
        model.viewObjectRepaint(this, output);
        return true;
    }

    public boolean processProfiles()
    {
        velocityModel = wizard.getVelocityModel();
        if (velocityModel == null) return true;
        panel.appendLine("Extracting profiles for the following model(s)...");
        panel.appendLine("\t" + velocityModel.getName());

        inputVelocityModel = wizard.getInputVelocityModel();
        if( inputVelocityModel != null)
            panel.appendLine("\t" + inputVelocityModel.getName());

        if(model == null) return false;

        StsPreStackLineSet currentLineSet = wizard.getCurrentLineSet();
        if (currentLineSet == null) return false;

        try
        {
            model.disableDisplay();

            StsCursor cursor = new StsCursor(model.getGlPanel3d(), Cursor.WAIT_CURSOR);
            logMessage("Preparing to extract profiles ....");

            if( velocityModel != null)
                success = doProcessProfile( velocityModel);

            if( success)
            {
                if( inputVelocityModel != null)
                    success = doProcessProfile( inputVelocityModel);
            } 

            panel.finished();
            panel.appendLine("Extraction Complete");
            panel.setDescriptionAndLevel("Extraction Complete", StsProgressBar.INFO);
            
            if( cursor!=null )
                cursor.restoreCursor();
            
            model.enableDisplay();
            
            if ( !success )
            {
                panel.setDescriptionAndLevel("Extraction Failed", StsProgressBar.FATAL);
            }
            else 
            {
                model.win3dDisplay();
            }

            enableFinish();
            return success;
       }
       catch (Exception e)
       {
           panel.appendLine("Extraction Failed");
           StsException.outputException("StsVelocityAnalysisWizard.extractVelocityProfiles() failed.", e, StsException.WARNING);
           panel.setDescriptionAndLevel("Exception thrown: " + e.getMessage(), StsProgressBar.ERROR);
           return false;
       }
    }

    public void run()
    {
        if (canceled)
        {
            success = false;
            return;
        }

        success = processProfiles();
    }

    public boolean end()
    {
        return true;
    }
}
