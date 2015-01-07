package com.Sts.Utilities.Interpolation;

import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jan 24, 2007
 * Time: 10:27:00 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsScatteredInterpolation
{
    int maxNObjects;
    StsScatteredDataObject[] weightedObjects = new StsScatteredDataObject[maxNObjects];
    int nObjects = 0;
    double maxDistSq = StsParameters.largeDouble;
    double weightSum = 0.0;

    public StsScatteredInterpolation(int maxNObjects)
    {
        this.maxNObjects = maxNObjects;
    }

    public boolean addObject(StsScatteredDataObject dataObject)
    {
        if(nObjects == maxNObjects)
        {
            int compare = dataObject.compareTo(weightedObjects[nObjects-1]);
            if(compare < 0)
                return replace(dataObject, nObjects-1);
            return false;
        }
        for (int n = nObjects - 1; n >= 0; n--)
        {
            int compare = dataObject.compareTo(weightedObjects[n]);
            if (compare >= 0)
            {
                if (nObjects >= maxNObjects - 2)
                    return false;
                else
                    return insert(dataObject, n + 1);
            }
        }
        return insert(dataObject, 0);
    }

    public boolean insert(StsScatteredDataObject object, int position)
    {

        if (position <= nObjects)
        {
            weightedObjects = (StsScatteredDataObject[]) StsMath.arrayInsertElementBefore(weightedObjects, object, position);
            nObjects++;
            return true;
        }
        return false;
    }

    public boolean replace(StsScatteredDataObject object, int position)
    {
        weightedObjects[position] = object;
        return true;
    }
    public boolean setRadialWeights()
    {
        double[] distances = new double[nObjects];
        for(int n = 0; n < nObjects; n++)
            distances[n] = weightedObjects[n].getDistance();
        double weightSum = 0.0;
        if(nObjects < 3)
        {
            for(int n = 0; n < nObjects; n++)
             {
                 double weight = 1.0 / (distances[n]);
                 weightSum += weight;
                 weightedObjects[n].weight = weight;
             }
        }
        else
        {
            double maxDistance = weightedObjects[nObjects-1].getDistance();
            for(int n = 0; n < nObjects; n++)
            {
                double dwt = (distances[n] - maxDistance);
                double weight = dwt * dwt / (distances[n] * maxDistance);
                weightSum += weight;
                weightedObjects[n].weight = weight;
            }
        }

        if(weightSum == 0.0)
        {
		    System.out.println("StsSpiralRadialInterpolation.Weights failed. weightSum = 0.0");
            return false;
        }
        for(int n = 0; n < nObjects; n++)
            weightedObjects[n].weight /= weightSum;
        return true;
    }

    public StsScatteredDataObject[] getWeightedObjects()
    {
        return (StsScatteredDataObject[])StsMath.trimArray(weightedObjects, nObjects);
    }

    public StsScatteredDataObject interpolate(StsScatteredDataObject result)
    {
        for(int n = 0; n < nObjects; n++)
            result = weightedObjects[n].interpolateValue(result);
        return result;
    }

    static public void main(String[] args)
    {
        StsScatteredInterpolation weighter = new StsScatteredInterpolation(3);
        WeightedObject[] weightedObjects = new WeightedObject[4];
        weightedObjects[0] = new WeightedObject(10, new Float(10));
        weightedObjects[1] = new WeightedObject(2, new Float(2));
        weightedObjects[2] = new WeightedObject(1, new Float(1));
        weightedObjects[3] = new WeightedObject(5, new Float(5));
        for(int n = 0; n < 4; n++)
            weighter.addObject(weightedObjects[n]);
        weighter.setRadialWeights();
        StsScatteredDataObject[] finalObjects = weighter.weightedObjects;
        WeightedObject result = new WeightedObject(0, new Float(0));
        result = (WeightedObject)weighter.interpolate(result);
/*        
        for(int n = 0; n < finalObjects.length; n++)
        {
            result = (WeightedObject)finalObjects[n].interpolateValue(result);
            float value = ((WeightedObject)finalObjects[n]).getFloatValue();
            double weight = finalObjects[n].weight;
            System.out.print(n + " value " + value + " weight " + weight + "   ");
        }
*/
        System.out.println("\n interpolatedValue " + result.getFloatValue());
    }
}

class WeightedObject extends StsScatteredDataObject
{
    WeightedObject(double distSq, Object value)
    {
        super(distSq, value);
    }

    public StsScatteredDataObject interpolateValue(StsScatteredDataObject result)
    {
        float value = getFloatValue();
        WeightedObject resultObject = (WeightedObject)result;
        float resultValue = resultObject.getFloatValue();
        resultObject.setFloatValue(resultValue + (float)weight*value);
        return resultObject;
      }

    public float getFloatValue()
    {
        return ((Float)valueObject).floatValue();
    }

    public Object setFloatValue(float value)
    {
        return new Float(value);
    }
}


