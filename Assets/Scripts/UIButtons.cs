using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class UIButtons : MonoBehaviour
{
    public GameObject manager; 
    public void StartExperience()
    {
        EventsPool.StartExperienceEvent.Invoke();
    }

    public void StartVolcanoMode()
    {
        manager.GetComponent<GPU_Particle_Manager>().showVolcanoIndicator = true; 
        
    }
    public void StartTsunamiMode()
    {
        manager.GetComponent<GPU_Particle_Manager>().showTsunamiIndicator = true;
    }
}
