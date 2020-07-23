from compass import models
from compass.models.MetabolicModel import MetabolicModel

import pytest


@pytest.mark.slow
def test_sbml_3():
    model = models.load_metabolic_model("RECON1_xml")
    assert isinstance(model, MetabolicModel)
    assert len(model.reactions) == 3742
    assert len(model.species) == 2766


@pytest.mark.slow
def test_sbml_2():
    model = models.load_metabolic_model("RECON2.2")
    assert isinstance(model, MetabolicModel)
    assert len(model.reactions) == 7785
    assert len(model.species) == 6047


@pytest.mark.slow
def test_mat():
    model = models.load_metabolic_model("RECON2_mat")
    assert isinstance(model, MetabolicModel)
    assert len(model.reactions) == 7440
    assert len(model.species) == 5063


@pytest.mark.slow
def test_to_json():
    model = models.load_metabolic_model("RECON2.2")
    json = model.to_JSON()
    assert isinstance(json, str)

    model = models.load_metabolic_model("RECON1_xml")
    json = model.to_JSON()
    assert isinstance(json, str)

    model = models.load_metabolic_model("RECON2_mat")
    json = model.to_JSON()
    assert isinstance(json, str)
