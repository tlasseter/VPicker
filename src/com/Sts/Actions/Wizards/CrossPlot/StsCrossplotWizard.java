package com.Sts.Actions.Wizards.CrossPlot;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

public class StsCrossplotWizard extends StsWizard
{
    private StsSelectCrossplot selectCrossplot;
    private StsDefineCrossplot defineCrossplot;

    private StsCrossplot selectedCrossplot = null;
    private StsTypeLibrary selectedLibrary = null;

    private boolean hasAttributeData = false;
    
    private StsWizardStep[] mySteps =
    {
        selectCrossplot = new StsSelectCrossplot(this),
        defineCrossplot = new StsDefineCrossplot(this)
    };

    public StsCrossplotWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Select/Define Crossplot");
        dialog.getContentPane().setSize(500, 600);
        if(!super.start()) return false;

        StsCrossplotClass crossplotClass = StsCrossplot.getCrossplotClass();
        int nCrossplots = crossplotClass.getSize();

        StsCrossplot currentCrossplot = crossplotClass.getCurrentCrossplot();
        if(currentCrossplot == null && nCrossplots > 0)
        {
            currentCrossplot = (StsCrossplot)crossplotClass.getLast();
            crossplotClass.setCurrentObject(currentCrossplot);
//            crossplotClass.setCurrentCrossplot(currentCrossplot);
        }
        selectedLibrary = StsTypeLibrary.getCreateGenericLibrary();

        return true;
    }

    public int getNumberOfAvailableVolumes()
    {
        int nSeismicVolumes = model.getStsClassSize(StsSeismicVolume.class);
        int nVirtualVolumes = ((StsVirtualVolumeClass)model.getCreateStsClass(StsVirtualVolume.class)).getVirtualVolumes().length;
        return nSeismicVolumes+nVirtualVolumes;    	
    }
    
    public boolean checkPrerequisites()
    {
        if(getNumberOfAvailableVolumes() < 2)
        {
            reasonForFailure = "Less than 2 seismic and/or virtual volumes loaded. Load additional volumes before continuing.";
            return false;
        }    	
    	return true;
    }
    
    public void checkAddToolbar()
    {
        if(currentStep == null)
            return;
        StsCrossplotToolbar crossplotToolbar = (StsCrossplotToolbar)model.win3d.getToolbarNamed(StsCrossplotToolbar.NAME);
        if(crossplotToolbar == null)
        {
            crossplotToolbar = new StsCrossplotToolbar(model.win3d);
            model.win3d.addToolbar(crossplotToolbar);
        }        
 //       StsViewSelectToolbar tb = model.win3d.getViewSelectToolbar();
 //       tb.setButtonVisibility(tb.getComponentNamed("XPView"), true);
    }

    public boolean end()
    {     	
        return super.end();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == selectCrossplot)
        {
            selectedCrossplot = selectCrossplot.getSelectedCrossplot();
            if(selectedCrossplot == null)
            {
                disableFinish();
                new StsMessage(frame, StsMessage.INFO, "No cross plot selected.");
                return;
            }
            else
            {
                enableFinish();
                model.setCurrentObjectDisplayAndToolbar(selectedCrossplot);
                model.setCurrentObject(selectedCrossplot);
                boolean changed = model.win3d.getCursor3d().setObject(selectedCrossplot);
//Redundant                StsComboBoxToolbar toolbar = (StsComboBoxToolbar)model.win3d.getToolbarNamed("Object Selection Toolbar");
//Redundant                toolbar.comboBoxSetItem(selectedCrossplot);
//                model.win3d.glPanel3d.checkAddView(StsViewXP.class);
//				model.setViewPreferred(StsViewXP.class, StsViewXP.viewClassnameXP);
//                model.win3d.getViewSelectToolbar().setButtonEnabled(model.glPanel3d.getCurrentView(), true);
//                 if(changed) model.win3d.glPanel3d.repaint();
                checkAddToolbar();
                hasAttributeData = selectedCrossplot.hasAttributeData();
                StsCrossplotToolbar crossplotToolbar = (StsCrossplotToolbar)model.win3d.getToolbarNamed(StsCrossplotToolbar.NAME);
                if(crossplotToolbar != null)
                	crossplotToolbar.enableDensityToggle(hasAttributeData);                                
                super.finish();
            }
        }
    }

    public void hasAttributeVolume(boolean val)
    {
    	hasAttributeData = val;
    }
    
    public void finish()
    {
        checkAddToolbar();       
        super.finish();
    }

    public StsSeismicVolume[] getSeismicVolumes()
    {
        StsSeismicVolume[] volumes = (StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class);
        StsVirtualVolume[] virtualVolumes = ((StsVirtualVolumeClass)model.getCreateStsClass(StsVirtualVolume.class)).getVirtualVolumes();

        if(virtualVolumes == null)
            return volumes;
        else if(virtualVolumes.length == 0)
            return volumes;
        else
            volumes = (StsSeismicVolume[])StsMath.arrayAddArray(volumes, virtualVolumes);
        return volumes;
    }

    public void setSelectedLibrary(StsTypeLibrary lib)
    {
        selectedLibrary = lib;
    }
    public StsTypeLibrary getSelectedLibrary(){ return selectedLibrary; }

    public void createNewCrossplot()
    {
        gotoStep(defineCrossplot);
        enableFinish();
    }

    static void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsCrossplotWizard crossplotWizard = new StsCrossplotWizard(actionManager);
        crossplotWizard.start();
    }
}
