using UnityEngine.Events;

public static class EventsPool
{
    public static UnityEvent StartExperienceEvent = new UnityEvent();

    public static UnityEvent ResetExperienceEvent = new UnityEvent();

    public static UnityEvent UpdateUIEvent = new UnityEvent();

}