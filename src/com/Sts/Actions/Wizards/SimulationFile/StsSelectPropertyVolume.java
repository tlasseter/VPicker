package com.Sts.Actions.Wizards.SimulationFile;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.RFUI.*;

import java.awt.*;

public class StsSelectPropertyVolume extends StsWizardStep
{
    StsSimulationModelPanel panel;
    StsHeaderPanel header;
    StsPropertyVolumeOld[] properties = null;

    public StsSelectPropertyVolume(StsWizard wizard)
    {
        super(wizard, new StsSimulationModelPanel(), null, new StsHeaderPanel());
        panel = (StsSimulationModelPanel) getContainer();
        panel.setPreferredSize(new Dimension(300,200));

        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Simulation File Export");
        header.setSubtitle("Select Property PostStack3d");
    }

    public boolean start()
    {
        // get the property volumes
    	if( properties == null )
        {
            properties = (StsPropertyVolumeOld[])model.getCastObjectList(StsPropertyVolumeOld.class);
            int nProperties = properties == null ? 0 : properties.length;
            if( nProperties > 0 )
            {
                String[] items = new String[nProperties];
                for( int i=0; i<nProperties; i++ )
                {
                    items[i] = new String(properties[i].getName());
                }
                panel.setPropertyModels(items);
            }
        }
        return properties != null;
    }


    public String trimName(String fullString)
    {
    	int index = fullString.indexOf("(");
        if( index > 0 )
			return fullString.substring(0, index-1);
        return fullString;
    }

    public boolean end() { return getSelectedPorosity() != null; }

    public StsPropertyVolumeOld getSelectedPorosity()
    {
        StsPropertyVolumeOld property = null;
		int item = panel.getPorosityModelIndex();
        if( item >= 0 ) return properties[item];
    	return null;
    }
    public StsPropertyVolumeOld getSelectedPermeability()
    {
        StsPropertyVolumeOld property = null;
		int item = panel.getPermeabilityModelIndex();
        if( item >= 0 ) return properties[item];
    	return null;
    }

    public double getKx() { return panel.getKx(); }
    public double getKy() { return panel.getKy(); }
    public double getKz() { return panel.getKz(); }
}
