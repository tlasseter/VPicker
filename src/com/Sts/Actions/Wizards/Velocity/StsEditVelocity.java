package com.Sts.Actions.Wizards.Velocity;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.StsProject;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

/**
 * <p>Title: jS2S development</p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author T.J.Lasseter
 * @version c51c
 */
public class StsEditVelocity extends StsWizardStep
{
	StsEditVelocityPanel panel;
	StsHeaderPanel header, info;

	/**
	 * StsEditVelocity creates or modifies the velocity model used for
	 * time-to-depth conversion.
	 *
	 * @param wizard StsVelocityWizard of which this is a step.
	 */
	public StsEditVelocity(StsWizard wizard)
	{
		super(wizard);
		panel = new StsEditVelocityPanel(this);
		header = new StsHeaderPanel();
		setPanels(panel, header);
//        int height = header.displayableHeight + panel.displayableHeight;
        panel.setPreferredSize(new Dimension(400, 400));
		header.setTitle("Velocity Model Construction");
		header.setSubtitle("Select Model Components");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Velocity");
		header.setInfoText(wizardDialog,"(1) Select the data components to use in the velocity model.\n"
                    + "  **** Optionally check box to use well marker-surface ties in model.\n"
                    + "  **** Optionally select one or more horizons from the list.\n"
                    + "  **** Optionally select imported velocity model from list.\n"
                    + "(2) Specify the depth datum.\n"
                    + "(3) Specify the time datum.\n"
                    + "(4) Press Next> Button\n"
                    + "  **** If only surfaces were selected, assign interval velocities is next.\n"
                    + "  **** If a volume was selected, volume definition is next.");
	}
	public boolean start()
	  {
          StsModelSurface[] modelSurfaces = (StsModelSurface[])model.getCastObjectList(StsModelSurface.class);
          StsSurface[] surfaces = StsSurface.getTimeSurfaces(modelSurfaces);
		  modelSurfaces = (StsModelSurface[])StsSurface.sortSurfaces(surfaces);
          panel.setSurfaces(modelSurfaces);
		  ((StsVelocityWizard)wizard).setSurfaces(modelSurfaces);

          StsSeismicVolume[] volumes = (StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class);
          panel.setVolumes(volumes);
          setTimeDepthDatum();

          int total = 0;
          if(volumes != null)
               total += volumes.length;
          if(surfaces != null)
               total += surfaces.length;
          if(total == 0)
              panel.initializePanelForSonics();

		  return true;
	  }

      private void setTimeDepthDatum()
      {
          float timeDatum = 0.0f;
          panel.topTimeDatumBean.setValue(timeDatum);
          ((StsVelocityWizard)wizard).setTopTimeDatum(timeDatum);
          StsObjectList wells = model.getInstances(StsWell.class);
          int nWells = wells.getSize();
          float minDepthDatum = StsParameters.largeFloat;
          float maxDepthDatum = -StsParameters.largeFloat;

          for(int n = 0; n < nWells; n++)
          {
              StsWell well = (StsWell)wells.getElement(n);
              StsLogCurve tdCurve = well.getLastLogCurveOfType(StsWellKeywordIO.TIME);
              if(tdCurve == null) continue;
              StsLogVector timeVector = tdCurve.getValueVector();
              float indexF = timeVector.getIndexF(0.0f);
              if(indexF == StsParameters.nullValue) continue;
              StsLogVector depthVector = tdCurve.getDepthVector();
              float depth = depthVector.getValue(indexF);
              minDepthDatum = Math.min(minDepthDatum, depth);
              maxDepthDatum = Math.max(maxDepthDatum, depth);
          }
          float depthDatum;
          if(minDepthDatum == StsParameters.largeFloat)
            depthDatum = 0.0f;
          else
          {
              if(maxDepthDatum - minDepthDatum > 1.0f)
              {
                  StsMessage.printMessage("Depth datum values computed from all td curves range from " + minDepthDatum + " to " + maxDepthDatum + ".\n" +
                                 "Set Top depth datum manually.");
              }
              depthDatum = (maxDepthDatum + minDepthDatum)/2;
          }
          panel.topDepthDatumBean.setValue(depthDatum);
          ((StsVelocityWizard)wizard).setTopDepthDatum(depthDatum);
     }

     public boolean end()
     {
         ((StsVelocityWizard)wizard).setSurfaces(panel.getSelectedSurfaces());
         return true;
     }
}
