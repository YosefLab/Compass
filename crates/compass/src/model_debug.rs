use gsmm::model::{self, Model};

pub fn list_reactions(model: Model) {
    for reaction in model.reactions() {
        println!("{}", reaction.name());
    }
}

pub fn reaction_info(model: Model, reaction: String) {
    let reaction = model
        .reactions()
        .iter()
        .find(|r| r.name() == reaction)
        .expect("Reaction not found");
    println!("{}", model.reaction_display(reaction));
}
