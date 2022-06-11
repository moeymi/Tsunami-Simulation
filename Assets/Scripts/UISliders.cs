using System.Collections;
using System.Collections.Generic;
using TMPro;
using UnityEngine;
using UnityEngine.UI;

public class UISliders : MonoBehaviour
{
    [SerializeField]
    GPU_Particle_Manager particleManager;

    [SerializeField]
    private Slider radiusSlider;

    [SerializeField]
    private Slider massSlider;

    [SerializeField]
    private Slider gcSlider;

    [SerializeField]
    private Slider restDensSlider;

    [SerializeField]
    private Slider viscositySlider;

    [SerializeField]
    private Slider dampingSlider;

    [SerializeField]
    private Slider particlesNumSlider;

    [SerializeField]
    private Slider dimsSlider;

    [SerializeField]
    private Slider volcanoInts;

    [SerializeField]
    private Slider volcanoRadius;

    [SerializeField]
    private Slider tsunamiInts;

    [SerializeField]
    private Slider tsunamiWidth;

    private TextMeshProUGUI radius;
    private TextMeshProUGUI mass;
    private TextMeshProUGUI gc;
    private TextMeshProUGUI resDens;
    private TextMeshProUGUI vc;
    private TextMeshProUGUI dmp;
    private TextMeshProUGUI prtNum;
    private TextMeshProUGUI dims;
    private TextMeshProUGUI volcanoIntsText;
    private TextMeshProUGUI volcanoRadiusText;
    private TextMeshProUGUI tsunamiIntsText;
    private TextMeshProUGUI tsunamiWidthText;
    private void Awake()
    {
        EventsPool.StartExperienceEvent.AddListener(InitializeValues);
        EventsPool.UpdateUIEvent.AddListener(InitializeValues);
        EventsPool.UpdateUIEvent.AddListener(UpdateUI);
        void updateUI(float f){
            EventsPool.UpdateUIEvent.Invoke();
        }
        try
        {
            radiusSlider.onValueChanged.AddListener(updateUI);
            massSlider.onValueChanged.AddListener(updateUI);
            gcSlider.onValueChanged.AddListener(updateUI);
            restDensSlider.onValueChanged.AddListener(updateUI);
            viscositySlider.onValueChanged.AddListener(updateUI);
            dampingSlider.onValueChanged.AddListener(updateUI);
            particlesNumSlider.onValueChanged.AddListener(updateUI);
        } catch { }
        try
        {
            volcanoInts.onValueChanged.AddListener(updateUI);
            volcanoRadius.onValueChanged.AddListener(updateUI);
            tsunamiInts.onValueChanged.AddListener(updateUI);
            tsunamiWidth.onValueChanged.AddListener(updateUI);
        }
        catch { }
    }
    private void UpdateUI()
    {
        if (radius == null)
        {
            try
            {
                radius = radiusSlider.transform.parent.Find("Value").GetComponent<TextMeshProUGUI>();
                mass = massSlider.transform.parent.Find("Value").GetComponent<TextMeshProUGUI>();
                gc = gcSlider.transform.parent.Find("Value").GetComponent<TextMeshProUGUI>();
                resDens = restDensSlider.transform.parent.Find("Value").GetComponent<TextMeshProUGUI>();
                vc = viscositySlider.transform.parent.Find("Value").GetComponent<TextMeshProUGUI>();
                dmp = dampingSlider.transform.parent.Find("Value").GetComponent<TextMeshProUGUI>();
                prtNum = particlesNumSlider.transform.parent.Find("Value").GetComponent<TextMeshProUGUI>();
                dims = dimsSlider.transform.parent.Find("Value").GetComponent<TextMeshProUGUI>();
            }
            catch { }
            try
            {
                volcanoIntsText = volcanoInts.transform.parent.Find("Value").GetComponent<TextMeshProUGUI>();
                volcanoRadiusText = volcanoRadius.transform.parent.Find("Value").GetComponent<TextMeshProUGUI>();
                tsunamiIntsText = tsunamiInts.transform.parent.Find("Value").GetComponent<TextMeshProUGUI>();
                tsunamiWidthText = tsunamiWidth.transform.parent.Find("Value").GetComponent<TextMeshProUGUI>();
            }
            catch { }
        }
        try
        {
            radius.text = radiusSlider.value.ToString();
            mass.text = massSlider.value.ToString();
            gc.text = gcSlider.value.ToString();
            resDens.text = restDensSlider.value.ToString();
            vc.text = viscositySlider.value.ToString();
            dmp.text = dampingSlider.value.ToString();
            prtNum.text = ((int)(particlesNumSlider.value - particlesNumSlider.value % 100)).ToString();
            dims.text = dimsSlider.value.ToString();
        } catch { }
        try
        {
            volcanoIntsText.text = volcanoInts.value.ToString();
            volcanoRadiusText.text = volcanoRadius.value.ToString();
            tsunamiIntsText.text = tsunamiInts.value.ToString();
            tsunamiWidthText.text = tsunamiWidth.value.ToString();
        }
        catch { }
    }

    private void InitializeValues()
    {
        try
        {
            particleManager.radius = radiusSlider.value;
            particleManager.mass = massSlider.value;
            particleManager.restDensity = restDensSlider.value;
            particleManager.viscosityCoefficient = viscositySlider.value;
            particleManager.damping = dampingSlider.value;
            int particlesNum = (int)(particlesNumSlider.value - particlesNumSlider.value % 100);
            particleManager.numberOfParticles = particlesNum;
            int y = 67108864;
            particleManager.maximumParticlesPerCell = Mathf.Min(200, (y / (particlesNum)));
            particleManager.dimensions = (int)dimsSlider.value;
        }
        catch { }
        try
        {
            particleManager.VolcanoRadius = volcanoRadius.value;
            particleManager.TsunamiWidth = tsunamiWidth.value;
        } catch { }
    }
}
